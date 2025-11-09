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

def get_type_by_name(name, resh_json):
    return list(
        filter(lambda x: x["name"] == name, resh_json["data-types"])
    )[0]

def test_ks():
    guid = "D4AAE7F3BE9448F6ED5F16DA1B97FFBD1"
    pdb = fetch_ms_pdb("ks.pdb", guid)
    resh_json = reshare_pdb.export(pdb).to_json_data()

    s_iosl = get_type_by_name("_IO_STACK_LOCATION",resh_json)
   
    assert s_iosl["size"] == 0x48
    print("[+] test_ks _IO_STACK_LOCATION size")

    assert len(s_iosl["content"]["members"]) == 9
    print("[+] test_ks _IO_STACK_LOCATION member count")

    s_owner_entry = get_type_by_name("_OWNER_ENTRY", resh_json)
    assert s_owner_entry["size"] == 0x10
    print("[+] test_ks _OWNER_ENTRY size")
    s_oe_union_name=s_owner_entry["content"]["members"][1]["type"]["type-name"]
    s_oe_union=get_type_by_name(s_oe_union_name, resh_json)
    assert len(s_oe_union["content"]["members"])==2
    print("[+] test_ks _OWNER_ENTRY union size")

test_ks()
