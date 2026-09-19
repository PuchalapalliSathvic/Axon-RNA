import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
import pandas as pd
import pytest
from axon_rna.real_data import (choose_transcript, validate_rna, sequence_features, parse_rnafold,
                               select_model_genes, annotation_features, encode_manual, assemble_dataset)
from axon_rna.real_modeling import make_estimator, make_splits


def test_transcript_preference_and_no_arbitrary_fallback():
    canonical={"id":"ENST1","MANE":[]}
    mane={"id":"ENST2","MANE":[{"type":"MANE_Select"}]}
    record={"canonical_transcript":"ENST1.5","Transcript":[canonical,mane]}
    assert choose_transcript(record)==(mane,"MANE_Select")
    record["Transcript"]=[canonical]
    assert choose_transcript(record)==(canonical,"Ensembl_canonical")
    record["canonical_transcript"]="ENST9"
    with pytest.raises(ValueError):
        choose_transcript(record)


def test_sequence_validation_does_not_delete_unknown_bases():
    assert validate_rna("aTgc")=="AUGC"
    assert sequence_features("AUGC")=={"sequence_length":4,"gc_content":.5}
    for seq in ["", "AUNG", "AU-G"]:
        with pytest.raises(ValueError):
            validate_rna(seq)


def test_fold_parser_and_normalization():
    result=parse_rnafold("GGGAAACCC\n(((...))) ( -3.60)\n","GGGAAACCC")
    assert result["normalized_mfe"]==pytest.approx(-.4)
    assert result["paired_fraction"]==pytest.approx(6/9)
    for output in ["((...))) (-2.0)", ").......( (-1.0)", "garbage"]:
        with pytest.raises(ValueError):
            parse_rnafold(output,"GGGAAACCC")


def test_dataset_selection_reproducible_and_ambiguous_positives_excluded():
    frame=pd.DataFrame({"transcript_id":["A","B","C","D"],"n_studies":[3]*4,
                        "direction_consistency":[1,1,2/3,2/3],"combined_p":[.001,.02,.5,.8],
                        "median_log2fc":[1,.2,.1,.1],"ad_dysregulated":[1,1,0,0]})
    candidates=pd.DataFrame({"transcript_id":["A"]})
    selected=select_model_genes(frame,candidates)
    pd.testing.assert_frame_equal(selected,select_model_genes(frame,candidates))
    assert len(selected)==2 and "B" not in selected.gene.tolist()
    assert selected.label.sum()==1
    with pytest.raises(ValueError):
        select_model_genes(pd.concat([frame,frame.iloc[:1]]),candidates)


def test_missing_review_is_not_negative_evidence():
    assert np.isnan(encode_manual("not reviewed"))
    assert np.isnan(encode_manual("unclear"))
    assert [encode_manual(v) for v in ["limited","moderate","strong"]]==[1,2,3]
    assert all(np.isnan(v) for v in annotation_features("","",False).values())
    assert annotation_features("synaptic signaling; RNA splicing","")["go_synaptic"]==1


def test_tfidf_fits_training_only_and_label_columns_are_rejected():
    frame=pd.DataFrame({"sequence":["AAAAAA","CCCCCC","AAAACC","CCCCAA"],"label":[0,1,0,1]})
    model=make_estimator([]).fit(frame,frame.label)
    vocabulary=model.named_steps["features"].named_transformers_["sequence"].vocabulary_
    assert "UUUU" not in vocabulary
    before=dict(vocabulary)
    model.predict_proba(pd.DataFrame({"sequence":["UUUUUU"]}))
    assert vocabulary==before
    with pytest.raises(ValueError):
        make_estimator(["evidence_score"])


def test_identical_sequences_stay_in_one_fold():
    rng=np.random.default_rng(7)
    seqs=["".join(rng.choice(list("ACGU"),size=40)) for _ in range(30)]
    frame=pd.DataFrame({"gene":[f"g{i}" for i in range(60)],"transcript_id":[f"t{i}" for i in range(60)],
                        "sequence":seqs*2,"label":([0,1]*15)*2})
    splits=make_splits(frame)
    held=[]
    for train,test in splits:
        assert not set(frame.sequence.iloc[train]) & set(frame.sequence.iloc[test])
        held.extend(test)
    assert sorted(held)==list(range(60))


def test_assembly_preserves_failures_and_rejects_many_to_many():
    selected=pd.DataFrame({"gene":["A","B"],"label":[1,0]})
    seq=pd.DataFrame({"gene":["A","B"],"transcript_id":["T1","T2"],"mapping_status":["mapped","failed"]})
    structure=pd.DataFrame({"gene":["A","B"],"transcript_id":["T1","T2"],"structure_status":["folded","mapping_failed"]})
    knowledge=pd.DataFrame({"gene":["A","B"],"go_count":[1,np.nan]})
    out=assemble_dataset(selected,seq,structure,knowledge)
    assert out.included_in_models.tolist()==[True,False]
    with pytest.raises(ValueError):
        assemble_dataset(selected,pd.concat([seq,seq]),structure,knowledge)


def test_scaler_and_imputer_are_fitted_only_to_training_values():
    frame=pd.DataFrame({"sequence":["AAAAAA","CCCCCC","AAAACC","CCCCAA"],
                        "normalized_mfe":[1.,2.,np.nan,3.],"label":[0,1,0,1]})
    model=make_estimator(["normalized_mfe"]).fit(frame,frame.label)
    numeric=model.named_steps["features"].named_transformers_["numeric"]
    assert numeric.named_steps["simpleimputer"].statistics_[0]==2
    assert numeric.named_steps["standardscaler"].mean_[0]==2
    model.predict_proba(pd.DataFrame({"sequence":["UUUUUU"],"normalized_mfe":[1000]}))
    assert numeric.named_steps["standardscaler"].mean_[0]==2
