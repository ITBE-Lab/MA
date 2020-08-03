
##
# @brief convert bytes to a NucSeq
# @details
# Usefull for converting reads stored as blob data in sqlite3 to NucSeq objects.
#
def nuc_seq_from_bytes(blob):
    converter = NucSeqSql()
    converter.fromBlob(blob)
    return converter.seq