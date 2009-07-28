typedef enum {Inheritance, Partof, Root} Strand_t;

typedef struct {
    int start,
    int end,
    Strand_t strand,
    int phase
} SeqLoc;

// Duplicates ID, Name, parents
typedef struct {
    char *tag,
    char **vals
} Attribute;

typedef struct {
    char *seqid,
    char *source,
    char *type,
    SeqLoc **locs,
    float score,
    char *ID,
    char *Name,
    char **parents,
    Attribute **attributes
} Feature;

typedef struct {
  char *ID,
  char *seqstr
} Seq;

typedef struct {
  Attribute **attributes,
  Feature **features,
  Seq **seqs
} GFFDoc;
