// Header file for gff_parser

// Author: Stephen W Veitch

// pragmas
#define GFF_FASTA 6
#define GFF_DOC 5
#define GFF_PRAGMA 4
// general success
#define GFF_SUCCESS 3
// comments
#define GFF_COMMENT 2
// eof
#define GFF_EOF 1
// fails
#define GFF_FAIL 0
#define GFF_NO_FILE -1
#define GFF_NO_MALLOC -2
#define GFF_BAD_RECORD -3
#define GFF_BAD_COLUMN -4
#define GFF_BAD_FASTA -5


typedef struct SeqLoc_list{
  SeqLoc *seqloc;
  struct SeqLoc_list *next;
} SeqLoc_list;

typedef struct Seq_list{
  Seq *seq;
  struct Seq_list *next;
} Seq_list;


int gff_get_doc(GFFDoc *gffdoc, char *file);
int gff_num_lines(char* file);
int gff_read_lines(int n, char *file);
int gff_features_match(Feature *f1, Feature *f2);
int gff_read_file(Feature **features, int *n, Seq **seqs, char *file);
int gff_read_fasta(Seq **seqs, FILE *fp);
int gff_read_pragma(FILE *fp);
int gff_read_record(Feature * feature, FILE *fp);
char **gff_split_attribute(char *attr, int *n);
//int gff_split_attribute(char **ret_values, char *attr, int *n);
char **gff_read_columns(char *record);
int gff_free_feature(Feature *f, int no_locs);
int gff_free(GFFDoc *gffdoc);
void gff_print_feature(Feature *f);
void gff_print(GFFDoc *gffdoc);


