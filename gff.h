typedef enum {Inheritance, Partof, Root} Strand_t;

typedef struct {
  int start;
  int end;
  //  Strand_t strand;
  int strand;
  int phase;
} SeqLoc;

// Duplicates ID, Name, parents
typedef struct {
  char *tag;
  char **vals;
  int num_vals;
} Attribute;

typedef struct {
  char *seqid;
  char *source;
  char *type;
  int num_locs;
  SeqLoc **locs;
  float score;
  char *ID;
  char *Name;
  int num_parents;
  char **parents;
  int num_attributes;
  Attribute **attributes;
} Feature;

typedef struct {
  char *ID;
  char *seqstr;
} Seq;

typedef struct {
  Attribute **attributes;
  int num_features;
  Feature **features;
  Seq **seqs;
} GFFDoc;
