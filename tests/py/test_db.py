from graph_id_db.finder import Finder


def test_hash():
    finder = Finder()
    docs = finder.find("7c22ad86c9eb7999")
    assert docs[0]["graph_id_type"] == "topo_md_id"
    assert docs[0]["proprietary_id"] == "mp-1183057"
    assert docs[0]["datasource"] == "MP"
    assert (
        docs[0]["url"]
        == "https://next-gen.materialsproject.org/materials/mp-1183057"
    )
    assert docs[0]["filehash"] == "ebefafe6539be2b6"
    assert len(docs) == 47


def test_not_found_hash():
    finder = Finder()
    docs = finder.find("xxxx")
    assert docs == []


def test_aflow():
    finder = Finder()

    fast_docs = finder.find("0000c690d6bd9179", is_fast=True)
    aflow_docs = finder.find("0000c690d6bd9179")

    assert fast_docs == []
    assert aflow_docs[0]["graph_id_type"] == "dc_id"
    assert aflow_docs[0]["proprietary_id"] == "aflow:06966bc7b36b2343"
    assert aflow_docs[0]["datasource"] == "AFLOW"
    assert (
        aflow_docs[0]["url"]
        == "https://aflow.org/material/?id=aflow:06966bc7b36b2343"
    )
    assert (
        aflow_docs[0]["filehash"]
        == "3ecb4cc4656f22e6c424f3fbaa0d2d4e62cfe1f827857072197ca3a31db16b78"
    )
    assert len(aflow_docs) == 2


def test_oqmd():
    finder = Finder()

    fast_docs = finder.find("00000278540685fb", is_fast=True)
    oqmd_docs = finder.find("00000278540685fb")

    assert fast_docs == []
    assert oqmd_docs[0]["graph_id_type"] == "dc_id"
    assert oqmd_docs[0]["proprietary_id"] == "573130"
    assert oqmd_docs[0]["datasource"] == "OQMD"
    assert oqmd_docs[0]["url"] == "https://oqmd.org/materials/entry/573130"
    assert (
        oqmd_docs[0]["filehash"]
        == "2da53b97cd761297533512ab3df4c400add0d6bf841585a1afff727f927ff171"
    )
    assert len(oqmd_docs) == 1


def test_pcod():
    finder = Finder()

    fast_docs = finder.find("00006e42109d6854", is_fast=True)
    pcod_docs = finder.find("00006e42109d6854")

    assert fast_docs == []
    assert pcod_docs[0]["graph_id_type"] == "md_id"
    assert pcod_docs[0]["proprietary_id"] == "8283166"
    assert pcod_docs[0]["datasource"] == "PCOD"
    assert (
        pcod_docs[0]["url"]
        == "https://www.crystallography.net/pcod/cif/8/828/8283166.cif"
    )
    assert pcod_docs[0]["filehash"] is None
    assert len(pcod_docs) == 2
