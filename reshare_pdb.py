# Based on GPL licensed code by Ry Auscitte
#
# https://auscitte.github.io/posts/Func-Prototypes-With-Pdbparse
# https://gist.github.com/Auscitte/37aa7b2d3be058cb6b4d5b8b4c13477a

from reshare import *
from reshare.helpers import *

import pdbparse
from pdbparse import tpi

from construct.core import *
import construct as cs

import json
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARNING)

# Parsers by Auscitte
ProcFrameData = cs.Struct(
    "rectyp"
    / cs.Enum(cs.Int16ul, S_FRAMEPROC=0x1012, S_CALLSITEINFO=0x1139, S_REGREL32=0x1111),
    "reminder"
    / cs.Switch(
        lambda ctx: ctx.rectyp,
        {
            "S_FRAMEPROC": "FRAMEPROCSYM"
            / cs.Struct(
                "cbFrame" / cs.Int32ul,
                "offPad" / cs.Int32ul,
                "cbSaveRegs" / cs.Int32ul,
                "offExHdlr" / cs.Int16ul,
                "flags" / cs.Int32ul,
            ),
            "S_REGREL32": "REGREL32"
            / cs.Struct(
                "off" / cs.Int32ul,
                "typind" / cs.Int32ul,
                "reg" / cs.Int16ul,
                "name" / cs.CString(encoding="utf8"),
            ),
            "S_CALLSITEINFO": "CALLSITEINFO"
            / cs.Struct(
                "off" / cs.Int32ul,
                "sect" / cs.Int16ul,
                "__reserved_0" / cs.Int16ul,
                "typind" / cs.Int32ul,
            ),
        },
    ),
)

ProcFrameEntries = cs.GreedyRange(
    cs.Struct(
        "reclen" / cs.Int16ul,
        "frame_entry"
        / cs.RestreamData(cs.Bytes(lambda ctx: ctx.reclen), ProcFrameData),
    )
)

GlobalProc = cs.Struct(
    "PROCSYM32"
    / cs.Struct(
        "reclen" / cs.Int16ul,
        "rectyp" / cs.Int16ul,
        "pParent" / cs.Int32ul,
        "pEnd" / cs.Int32ul,
        "pNext" / cs.Int32ul,
        "len" / cs.Int32ul,
        "DbgStart" / cs.Int32ul,
        "DbgEnd" / cs.Int32ul,
        "typind" / cs.Int32ul,
        "offset" / cs.Int32ul,
        "seg" / cs.Int16ul,
        "flags" / cs.Int8ul,
        "name" / cs.CString(encoding="utf8"),
    ),
    # making sure the entire length of PROCSYM32 has been parsed
    cs.Padding(lambda ctx: ctx.PROCSYM32.reclen + cs.Int16ul.sizeof() - ctx._io.tell()),
    "frame_data"
    / cs.RestreamData(
        # ctx.PROCSYM32.pEnd points to the region immediately following the last element of ProcFrameEntries
        # ctx.PROCSYM32.reclen does not include the reclen field hence the cs.Int16ul.sizeof() correction
        cs.Bytes(
            lambda ctx: ctx.PROCSYM32.pEnd
            - ctx._params.entry_offest
            - ctx.PROCSYM32.reclen
            - cs.Int16ul.sizeof()
        ),
        ProcFrameEntries,
    ),
)
GlobalProcSym = "PROCSYM32" / cs.Struct(
    "reclen" / cs.Int16ul,
    "rectyp" / cs.Int16ul,
    "pParent" / cs.Int32ul,
    "pEnd" / cs.Int32ul,
    "pNext" / cs.Int32ul,
    "len" / cs.Int32ul,
    "DbgStart" / cs.Int32ul,
    "DbgEnd" / cs.Int32ul,
    "typind" / cs.Int32ul,
    "offset" / cs.Int32ul,
    "seg" / cs.Int16ul,
    "flags" / cs.Int8ul,
    "name" / cs.CString(encoding="utf8"),
)
S_PROCREF = 0x1125


# base_type_size has been borrowed from pdbparse's pdb_print_ctypes.py
# (see https://github.com/moyix/pdbparse)
base_type_size = {
    "T_32PRCHAR": 4,
    "T_32PUCHAR": 4,
    "T_32PULONG": 4,
    "T_32PUQUAD": 4,
    "T_32PUSHORT": 4,
    "T_32PVOID": 4,
    "T_32PLONG": 4,
    "T_64PRCHAR": 8,
    "T_64PUCHAR": 8,
    "T_64PULONG": 8,
    "T_64PUQUAD": 8,
    "T_64PUSHORT": 8,
    "T_64PVOID": 8,
    "T_64PLONG": 8,
    "T_64PWCHAR": 8,
    "T_64PQUAD": 8,
    "T_64PSHORT": 8,
    "T_64PUINT": 8,
    "T_64PUINT4": 8,
    "T_64PINT4": 8,
    "T_64PCHAR": 8,
    "T_INT4": 4,
    "T_INT8": 8,
    "T_LONG": 4,
    "T_QUAD": 8,
    "T_RCHAR": 1,
    "T_REAL32": 4,
    "T_REAL64": 8,
    "T_REAL80": 10,
    "T_SHORT": 2,
    "T_UCHAR": 1,
    "T_UINT4": 4,
    "T_ULONG": 4,
    "T_UQUAD": 8,
    "T_USHORT": 2,
    "T_WCHAR": 2,
    "T_BOOL08": 8,
    "T_VOID": 0,
    "T_NOTYPE": 0,
    "T_HRESULT": 4,
}


types = {
    "ptr": ReshDataTypePy(
        name="void *",
        content=ReshDataTypeContentPrimitivePy(),
        size=8,
    ),
    "void *": ReshDataTypePy(
        name="void *",
        content=ReshDataTypeContentPrimitivePy(),
        size=8,
    ),
}


resh_union_count = 0


def is_overlapping(map, offset):
    for k, v in map.items():
        if offset >= k and offset < k + v[1]:
            return k
    return None

def create_structure(T):
    global resh_union_count

    ret = ReshDataTypePy(name=T.name, content=None, size=T.size)
    types[T.tpi_idx] = ret  # We cache the struct first to handle recursive structures

    content = ReshDataTypeContentStructurePy(members=[])
    if T.size == 0:
        ret.content = content
    else:
        last_bitfield = None
        struct_map = (
            {}
        )  # https://www.vergiliusproject.com/kernels/x64/windows-11/24h2/_DISPATCHER_HEADER
        substructs = list(filter(
            lambda x: hasattr(x, "index") and hasattr(x, "offset"),
            T.fieldlist.substructs,
        ))
        substructs.sort(key=lambda x: x.offset)
        for member in substructs:
            try:
                resh_member = None
                member_name = member.name
                if (
                    hasattr(member.index, "leaf_type")
                    and member.index.leaf_type == "LF_BITFIELD"
                ):
                    if last_bitfield != int(member.index.leaf_type):
                        resh_member = get_single_type(member.index.base_type)
                        last_bitfield = int(member.index.leaf_type)
                        member_name = "bitfield%X" % (last_bitfield)
                    else:
                        continue
                else:
                    last_bitfield = None
                    resh_member = get_single_type(member.index)

                overlapping_offset = is_overlapping(struct_map, member.offset)
                if overlapping_offset is not None:
                    resh_member_wrapped = ReshStructureMemberPy(
                        type=resh_member.name,
                        name=member_name,
                        offset=member.offset,
                    )

                    existing, _ = struct_map[overlapping_offset]
                    if existing.name.startswith("resh_union"):
                        types[existing.name].content.members.append(resh_member_wrapped)
                        if resh_member.size > struct_map[overlapping_offset][1]:
                            struct_map[overlapping_offset][1] = resh_member.size
                    else:
                        union_name = "resh_union%04X" % (resh_union_count)
                        resh_union_count += 1

                        size = struct_map[overlapping_offset][1]
                        if resh_member.size > size:
                            size = resh_member.size

                        union_content = ReshDataTypeContentUnionPy(
                            members=[existing, resh_member_wrapped]
                        )
                        union_type = ReshDataTypePy(
                            name=union_name,
                            content=union_content,
                            size=size,
                        )
                        types[union_name] = union_type
                        struct_map[member.offset] = [
                            ReshStructureMemberPy(
                                type=union_name,
                                name=union_name,
                                offset=member.offset,
                            ),
                            size,
                        ]
                else:
                    struct_map[member.offset] = [
                        ReshStructureMemberPy(
                            type=resh_member.name,
                            name=member_name,
                            offset=member.offset,
                        ),
                        resh_member.size,
                    ]
            except AttributeError:
                # print(dir(member))
                # print(member.name)
                # print(resh_member)
                logger.warning("Fucky member in struct: %s" % (T.name))
                raise
                continue

        content.members = [x[0] for _, x in sorted(struct_map.items())]
        ret.content = content
    return ret


def create_union(T):
    ret = ReshDataTypePy(name=T.name, content=None, size=T.size)
    content = ReshDataTypeContentUnionPy(members=[])
    try:
        if hasattr(T.fieldlist, "substructs"):
            last_bitfield = None
            for member in T.fieldlist.substructs:
                if not hasattr(member, "index"):
                    continue

                resh_member = None
                member_name = member.name
                if (
                    hasattr(member.index, "leaf_type")
                    and member.index.leaf_type == "LF_BITFIELD"
                ):
                    # **TODO** leaf_type is the same for all bitfields
                    # so this is not the correct way to distinguish bitfields
                    # that follow each other.
                    # Bitfields don't necessary include as many bits as their
                    # base type either
                    if last_bitfield != int(member.index.leaf_type):
                        resh_member = get_single_type(member.index.base_type)
                        last_bitfield = int(member.index.leaf_type)
                        member_name = "bitfield%X" % (last_bitfield)
                    else:
                        continue
                else:
                    last_bitfield = None
                    resh_member = get_single_type(member.index)

                size = resh_member.size

                offset = None
                if hasattr(member, "offset"):
                    offset = member.offset
                else:
                    offset = -1
                content.members.append(
                    ReshStructureMemberPy(
                        type=resh_member.name,
                        name=member_name,
                        offset=offset,
                    )
                )
        ret.content = content
        types[T.tpi_idx] = ret
    except AttributeError as ae:
        logger.warning("Empty union: %s" % (ret.name))
        raise
    return ret


def get_single_type(T):
    ret = None

    if "tpi_idx" not in dir(T):
        if str(T) in types:
            return types[str(T)]
        ret = ReshDataTypePy(
            name=str(T),
            content=ReshDataTypeContentPrimitivePy(),
            size=base_type_size[str(T)],
        )
        types[str(T)] = ret
    elif T.tpi_idx in types:
        return types[T.tpi_idx]
    elif T.leaf_type == "LF_BITFIELD":
        return None
    elif T.leaf_type == "LF_POINTER":
        pointed = get_single_type(T.utype)
        if pointed is None:
            raise Exception("Can't point to type:" + T.utype)
        content = ReshDataTypeContentPointerPy(
            target_type=pointed.name
        )
        ret = ReshDataTypePy(
            content=content, name=pointed.name + " *", size=8
        )
        types[T.tpi_idx] = ret
    elif T.leaf_type == "LF_STRUCTURE":
        ret = create_structure(T)
    elif T.leaf_type == "LF_UNION":
        ret = create_union(T)
    elif T.leaf_type == "LF_ENUM":
        # content=ReshDataTypeContentEnum()
        members = []
        base = get_single_type(T.utype)
        for f in T.fieldlist.substructs:
            if hasattr(f, "name"):
                name = f.name
                value = f.enum_value
                members.append(ReshEnumMember(name=name, value=value))
        ret = ReshDataTypePy(
            name=T.name,
            content=ReshDataTypeContentEnumPy(
                members=members,
                base_type=base.name,
            ),
            size=base.size,
        )
        types[T.tpi_idx] = ret
    elif T.leaf_type == "LF_MODIFIER":
        ret = get_single_type(T.modified_type)
        s = [mod for mod in ["const", "volatile", "unaligned"] if T.modifier[mod]]
        ret.modifiers = s
        types[T.tpi_idx] = ret
    elif T.leaf_type == "LF_ARRAY":
        member_type = get_single_type(T.element_type)
        member_type_size = member_type.size

        if member_type_size is None:
            print(member_type)
            raise Exception("Can't determine size for '%s'" % (member_type.name))

        if T.size % member_type_size != 0:
            raise Exception(
                "Not dividers: '%s' (%d) / '%s' (%d) "
                % (T.name, T.size, member_type.name, member_type_size)
            )

        array_item_count = int(T.size / member_type_size)
        content = ReshDataTypeContentArrayPy(
            base_type=ReshTypeSpec(type_name=member_type.name, embedded_type=member_type),
            length=array_item_count,
        )
        ret = ReshDataTypePy(
            content=content,
            name="%s[%d]" % (member_type.name, array_item_count),
            size=T.size,
        )
        types[T.tpi_idx] = ret
    elif T.leaf_type == "LF_NESTTYPE":
        ret = get_single_type(T.index)
        types[T.tpi_idx] = ret
    elif T.leaf_type == "LF_PROCEDURE":
        ret_type = get_single_type(T.return_type)
        arg_list = []
        if T.parm_count > 0:
            for i, arg in enumerate(T.arglist.arg_type):
                my_arg = get_single_type(arg)
                arg_list.append(
                    ReshFunctionArgumentPy(name="param%d" % (i,), type=my_arg.name)
                )

        content = ReshDataTypeContentFunctionPy(
            arguments=arg_list,
            return_type=ret_type.name,
            calling_convention=None,
        )
        ret = ReshDataTypePy(
            content=content,
            name="function%X" % (T.tpi_idx),
            size=0,
        )
        types[T.tpi_idx] = ret
        return ret
    elif T.leaf_type == "LF_CLASS" or T.leaf_type == "LF_VTSHAPE":
        return types["ptr"]

    return ret


def export_types(pdb):
    for idx, T in pdb.STREAM_TPI.types.items():
        if idx not in types:
            get_single_type(T)


def address_to_resh(offset):
    return ReshAddress(
        list(offset.to_bytes(8, byteorder="little", signed=False)),
        "ram",
    )


def export_function_symbols(pdb):
    f_map = {}
    for f in pdb.STREAM_GSYM.globals:
        if hasattr(f, "name") and hasattr(f, "offset"):
            f_map[f.name] = ReshSymbolPy(
                name=f.name,
                address=address_to_resh(f.offset),
                confidence=ReshSymbolConfidence.FACT,
            )

    fncs = list(filter(lambda s: s.leaf_type == S_PROCREF, pdb.STREAM_GSYM.globals))
    ret = []
    for f in fncs:
        data = pdb.streams[pdb.STREAM_DBI.DBIExHeaders[f.iMod - 1].stream].data
        x = GlobalProcSym.parse(data[f.offset :], entry_offest=f.offset)
        f_map[f.name].type = ReshTypeSpec(type_name="function%X" % (x.typind,), embedded_type=None)
    return f_map.values()


def export(pdb):
    export_types(pdb)
    f_symbols = export_function_symbols(pdb)

    resh = ResharePy(pdb.fp.name)

    for _, dt in types.items():
        resh.data_types.append(dt)

    resh.symbols.extend(f_symbols)
    return resh


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="The PDB file")

    args = parser.parse_args()
    pdb = pdbparse.parse(args.pdb)

    resh = export(pdb)

    print(json.dumps(resh.to_json_data(), indent=2))
