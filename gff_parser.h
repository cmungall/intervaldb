// Header file for gff_parser

GFFDoc *gff_get_doc(char *file);
int gff_num_lines(char* file);
int gff_read_lines(int n,char *file);
Feature **gff_read_file(int r,char *file);
Feature *gff_read_record(FILE *fp);
char **gff_split_attribute(char *attr, int *n);
char **gff_read_columns(char *record);
void gff_print_feature(Feature *f);


