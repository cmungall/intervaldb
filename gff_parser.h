// Header file for gff_parser

// Author: Stephen W Veitch

// pragmas
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

int gff_get_doc(GFFDoc *gffdoc, char *file);
int gff_num_lines(char* file);
int gff_read_lines(int n, char *file);
int gff_read_file(Feature **features, int *n, char *file);
int gff_read_pragma(FILE *fp);
int gff_read_record(Feature * feature, FILE *fp);
char **gff_split_attribute(char *attr, int *n);
//int gff_split_attribute(char **ret_values, char *attr, int *n);
char **gff_read_columns(char *record);
int gff_free(GFFDoc *gffdoc);
void gff_print_feature(Feature *f);


