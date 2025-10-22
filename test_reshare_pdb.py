import httpx
import reshare_pdb
import typing
import os
import pdbparse


samples_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "samples")


def fetch_ms_pdb(name: str, guid: str) -> pdbparse.PDBStream:
    pdb_path = os.path.join(samples_path, name)
    if not os.path.exists(pdb_path):
        content = httpx.get(
            "http://msdl.microsoft.com/download/symbols/%s/%s/%s" % (name, guid, name),
            follow_redirects=True,
        )
        print(content.status_code)
        with open(pdb_path, "wb") as out:
            out.write(content.content)

    return pdbparse.parse(pdb_path)


def test_ks():
    guid = "D4AAE7F3BE9448F6ED5F16DA1B97FFBD1"
    pdb = fetch_ms_pdb("ks.pdb", guid)
    resh_json = reshare_pdb.export(pdb).to_json_data()
    s_iosl = list(
        filter(lambda x: x["name"] == "_IO_STACK_LOCATION", resh_json["data-types"])
    )[0]

    assert s_iosl["size"] == 0x48
    print("[+] test_ks _IO_STACK_LOCATION size")

    assert len(s_iosl["content"]["members"]) == 9
    print("[+] test_ks _IO_STACK_LOCATION member count")


test_ks()
