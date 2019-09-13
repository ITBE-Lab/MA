#ifndef MINIMAP2_H
#define MINIMAP2_H

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define MM_F_NO_DIAG       0x001 // no exact diagonal hit
#define MM_F_NO_DUAL       0x002 // skip pairs where query name is lexicographically larger than target name
#define MM_F_CIGAR         0x004
#define MM_F_OUT_SAM       0x008
#define MM_F_NO_QUAL       0x010
#define MM_F_OUT_CG        0x020
#define MM_F_OUT_CS        0x040
#define MM_F_SPLICE        0x080 // splice mode
#define MM_F_SPLICE_FOR    0x100 // match GT-AG
#define MM_F_SPLICE_REV    0x200 // match CT-AC, the reverse complement of GT-AG
#define MM_F_NO_LJOIN      0x400
#define MM_F_OUT_CS_LONG   0x800
#define MM_F_SR            0x1000
#define MM_F_FRAG_MODE     0x2000
#define MM_F_NO_PRINT_2ND  0x4000
#define MM_F_2_IO_THREADS  0x8000
#define MM_F_LONG_CIGAR    0x10000
#define MM_F_INDEPEND_SEG  0x20000
#define MM_F_SPLICE_FLANK  0x40000
#define MM_F_SOFTCLIP      0x80000
#define MM_F_FOR_ONLY      0x100000
#define MM_F_REV_ONLY      0x200000
#define MM_F_HEAP_SORT     0x400000
#define MM_F_ALL_CHAINS    0x800000
#define MM_F_OUT_MD        0x1000000
#define MM_F_COPY_COMMENT  0x2000000
#define MM_F_EQX           0x4000000 // use =/X instead of M
#define MM_F_PAF_NO_HIT    0x8000000 // output unmapped reads to PAF
#define MM_F_NO_END_FLT    0x10000000

#define MM_I_HPC          0x1
#define MM_I_NO_SEQ       0x2
#define MM_I_NO_NAME      0x4

#define MM_IDX_MAGIC   "MMI\2"

#define MM_MAX_SEG       255

#ifdef __cplusplus
extern "C" {
#endif

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

// minimap2 index
typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
} mm_idx_seq_t;

typedef struct {
	int32_t b, w, k, flag;
	uint32_t n_seq;            // number of reference sequences
	int32_t index;
	mm_idx_seq_t *seq;         // sequence name, length and offset
	uint32_t *S;               // 4-bit packed sequence
	struct mm_idx_bucket_s *B; // index (hidden)
	void *km, *h;
} mm_idx_t;

// indexing and mapping options
typedef struct {
	short k, w, flag, bucket_bits;
	int mini_batch_size;
	uint64_t batch_size;
} mm_idxopt_t;

// index reader
typedef struct {
	int is_idx, n_parts;
	int64_t idx_size;
	mm_idxopt_t opt;
	FILE *fp_out;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;
} mm_idx_reader_t;

// global variables
extern int mm_verbose, mm_dbg_flag; // verbose level: 0 for no info, 1 for error, 2 for warning, 3 for message (default); debugging flag
extern double mm_realtime0; // wall-clock timer


/**
 * Initialize an index reader
 *
 * @param fn         index or fasta/fastq file name (this function tests the file type)
 * @param opt        indexing parameters
 * @param fn_out     if not NULL, write built index to this file
 *
 * @return an index reader on success; NULL if fail to open _fn_
 */
mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out);

/**
 * Read/build an index
 *
 * If the input file is an index file, this function reads one part of the
 * index and returns. If the input file is a sequence file (fasta or fastq),
 * this function constructs the index for about mm_idxopt_t::batch_size bases.
 * Importantly, for a huge collection of sequences, this function may only
 * return an index for part of sequences. It needs to be repeatedly called
 * to traverse the entire index/sequence file.
 *
 * @param r          index reader
 * @param n_threads  number of threads for constructing index
 *
 * @return an index on success; NULL if reaching the end of the input file
 */
mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads);

/**
 * Destroy/deallocate an index reader
 *
 * @param r          index reader
 */
void mm_idx_reader_close(mm_idx_reader_t *r);

int mm_idx_reader_eof(const mm_idx_reader_t *r);

/**
 * Check whether the file contains a minimap2 index
 *
 * @param fn         file name
 *
 * @return the file size if fn is an index file; 0 if fn is not.
 */
int64_t mm_idx_is_idx(const char *fn);

/**
 * Load a part of an index
 *
 * Given a uni-part index, this function loads the entire index into memory.
 * Given a multi-part index, it loads one part only and places the file pointer
 * at the end of that part.
 *
 * @param fp         pointer to FILE object
 *
 * @return minimap2 index read from fp
 */
mm_idx_t *mm_idx_load(FILE *fp);

/**
 * Append an index (or one part of a full index) to file
 *
 * @param fp         pointer to FILE object
 * @param mi         minimap2 index
 */
void mm_idx_dump(FILE *fp, const mm_idx_t *mi);

/**
 * Create an index from strings in memory
 *
 * @param w            minimizer window size
 * @param k            minimizer k-mer size
 * @param is_hpc       use HPC k-mer if true
 * @param bucket_bits  number of bits for the first level of the hash table
 * @param n            number of sequences
 * @param seq          sequences in A/C/G/T
 * @param name         sequence names; could be NULL
 *
 * @return minimap2 index
 */
mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name);

/**
 * Print index statistics to stderr
 *
 * @param mi         minimap2 index
 */
void mm_idx_stat(const mm_idx_t *idx);

/**
 * Destroy/deallocate an index
 *
 * @param r          minimap2 index
 */
void mm_idx_destroy(mm_idx_t *mi);


#ifdef __cplusplus
}
#endif

#endif // MINIMAP2_H
