// A parser for GFF3
//
// Author: Stephen W Veitch

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <limits.h>

#include "gff.h"
#include "gff_parser.h"

#define BUFFER_SIZE 1024

main()
{
  int err;
  int x;
  char file[BUFFER_SIZE];

  //  for (x=0; x<2; x++) 
  for (x=1; x<2; x++) 
    {
      x = 7;
      sprintf(file,"t/gfftab%d.dat",x);
      GFFDoc *gff = malloc(sizeof(GFFDoc));

      //err = gff_get_doc(gff,"t/gfftab.dat"); 
      err = gff_get_doc(gff,file); 

      if (err) // positive value so ok
	{
	  gff_print(gff);
	  //printf("Print file\n");
	  //for (x=0; x<gff->num_features; x++)
	  //gff_print_feature(gff->features[x]);
	}
      err = gff_free(gff);
      printf("X=%d\n",x);
    }
}

int gff_get_doc(GFFDoc *gffdoc, char *file)
{
  int x = 0;
  int r;
  int err;

  // count the number of lines and create an array of Features
  int n = gff_num_lines(file);
  printf("number of lines = %d\n",n);
  
  //r = gff_read_lines(n,file);
  //printf("number of records = %d\n",r);

  // TODO sort out below to get right amount of memory allocated

  // allocate enough memory for more records than we need
  gffdoc->features = malloc(sizeof(Feature*)*n);
  gffdoc->seqs = malloc(sizeof(Seq*)*n);
  err = gff_read_file(gffdoc->features, &gffdoc->num_features, gffdoc->seqs, &gffdoc->num_seqs, file);

  // handle errors
  switch(err)
    {
    case GFF_FAIL:
      printf("ERROR: General failure\n");
      break;
    case GFF_NO_FILE:
      printf("ERROR: Failure to open file %s\n",file);
      break;
    case GFF_NO_MALLOC:
      printf("ERROR: Failure to allocate memory\n");
      break;
    }

  return err;
}

// reads the number of lines in the file
// NB should be refined as the number of features allocates
// is dependant on this but not necessary
int gff_num_lines(char* file)
{
  FILE *fp;
  int n = 0;
  int c;

  // open file and count lines
  if (!(fp = fopen(file,"r"))) return GFF_NO_FILE;

  while ((c = fgetc(fp)) != EOF)
    {
      if (c == '\n') n++;
      //      printf("%c",c);  
    }

  fclose(fp);
  return n;
}

int gff_read_lines(int n,char *file)
{
  FILE *fp;
  int i = 0;
  int r = 0;
  char line[100];

  // open file and count lines
  if (!(fp = fopen(file,"r"))) return GFF_NO_FILE;

  while (i<n && fscanf(fp,"%s\n",&line))
    {
      i++;
      // check line is not a comment or empty
      if (line[0] != '#' && line[0] != '\n') 
	{
	r++;
	//printf("%s\n",line);
	}
    }

  fclose(fp);
  return r;
}

// returns positive if features match
int gff_features_match(Feature *f1, Feature *f2)
{
  printf("features match:1\n");

  printf("features match <%s><%s>\n",f1->seqid,f2->seqid);
  if (!(strcmp(f1->seqid,f2->seqid)==0 &&
	strcmp(f1->source,f2->source)==0 &&
	strcmp(f1->type,f2->type)==0 &&
	strcmp(f1->ID,f2->ID)==0 &&
	strcmp(f1->Name,f2->Name)==0))
    return 0;

  int n,v;
  for (n=0; n<f1->num_parents; n++)
    if (!(strcmp(f1->parents[n],f2->parents[n])==0))
	return 0;
  for (n=0; n<f1->num_attributes; n++)
    if ((strcmp(f1->attributes[n]->tag,f2->attributes[n]->tag)==0) &&
      f1->attributes[n]->num_vals==f1->attributes[n]->num_vals)
      {
	for (v=0; v<f1->attributes[n]->num_vals; v++)
	  if (!(strcmp(f1->attributes[n]->vals[v],f2->attributes[n]->vals[v])))
	    return 0;
      }	
    else
      return 0;
  
  printf("features matched\n");
  return 1;
}

int gff_read_file(Feature **features,int *n, Seq **seqs, int *num_seqs, char *file)
{
  FILE *fp;
  int x = 0;
  Feature *feature;
  int err = GFF_SUCCESS;
  SeqLoc_list *seqloc_list;
  SeqLoc_list *seqloc_p;
  int matched = 0;
  int match = 0;

  // open file and count lines
  if (!(fp = fopen(file,"r"))) return GFF_NO_FILE;
  
  // first check that we are dealing with a GFF3 file
  feature = malloc(sizeof(Feature)); 
  if ((err = gff_read_record(feature,fp)) != GFF_DOC)
    {
      free(feature);
      return GFF_FAIL;
    }

  while (err != GFF_EOF && err != GFF_FASTA)
    {
      printf("gff_read_file:1 x=%d\n",x);
      feature = malloc(sizeof(Feature));
      feature->num_locs = 1;

      err = gff_read_record(feature,fp);
      match = 0;

      switch(err)
	{
	case GFF_SUCCESS: 
	  printf("gff_read_file:1 x=%d\n",x);
	  features[x++] = feature;
	  // check if the feature matches with previous one
	  if (x>1 && gff_features_match(features[x-2],feature))
	    {
	      match = 1;
	      printf("gff_read_file:1a x=%d start=%d\n",x,*feature->locs[0]->start);
	      if (matched == 0)
		{
		  seqloc_list = (SeqLoc_list*)malloc(sizeof(SeqLoc_list));
		  seqloc_list->seqloc = features[x-2]->locs[0];
		  seqloc_p = seqloc_list;
		  matched++;
		}
	      seqloc_p->next = (SeqLoc_list*)malloc(sizeof(SeqLoc_list));
	      seqloc_p = seqloc_p->next;
	      seqloc_p->seqloc = feature->locs[0];
	      seqloc_p->next = NULL;
	      matched++;
	      //printf("gff_read_file:1b x=%d\n",x);
	      gff_free_feature(feature,1);
	      x--;
	      printf("gff_read_file:1c x=%d matched=%d\n",x,matched);
	    }
	  break;
	case GFF_PRAGMA: 
	  printf("line %d is a PRAGMA\n",x);
	  free(feature);
	  break;
	case GFF_COMMENT: 
	  printf("line %d is a Comment\n",x);
	  free(feature);
	  break;
	case GFF_FASTA: 
	  printf("Remainder is FASTA data\n",x);
	  // must deal with this but currently just stop
	  err = gff_read_fasta(seqs,num_seqs,fp);
	  if (err == GFF_BAD_FASTA)
	    {
	      free(feature);
	      return GFF_BAD_FASTA;
	    }
	  //err = GFF_EOF;
	  //err = GFF_FASTA;
	  //free(feature);
	  break;
	}

      printf("matched=%d match=%d\n",matched,match);
      // we have the end of a matched features
      if (matched && !match)
	{ // add the list of SeqLoc's to the feature
	  if (err == GFF_EOF) x++; //this is for a last matched feature!
	  match = 0;
	  printf("gff_read_file:1d x=%d matched=%d\n",x,matched);
	  free(features[x-2]->locs);
	  features[x-2]->locs = malloc(sizeof(SeqLoc*)*matched);
	  features[x-2]->num_locs = matched;
	  printf("gff_read_file:1e x=%d matched=%d\n",x,matched);
	    
	  seqloc_p = seqloc_list;
	  printf("start=%d num_locs=%d\n",*(seqloc_p->seqloc->start),features[x-2]->num_locs);
	  int m;
	  printf("gff_read_file:1f x=%d matched=%d\n",x,matched);
	  for (m=0; m<matched; m++)
	    {
	      printf("gff_read_file:1g x=%d match=%d\n",x,m);
	      features[x-2]->locs[m] = seqloc_p->seqloc;
	      printf("gff_read_file:1h x=%d match=%d\n",x,m);
	      printf("matched %d start=%d\n",m,*(seqloc_p->seqloc->start));
	      seqloc_p = seqloc_list->next;
	      free(seqloc_list);
	      seqloc_list = seqloc_p;
	    }
	  matched = 0;
	  if (err == GFF_EOF) x--; //this is for a last matched feature!
	}
    }

  free(feature);
  *n = x;

  printf("gff_read_file:2\n");
  fclose(fp);
  printf("gff_read_file:3\n");

  return err;
}

int gff_read_fasta(Seq **seqs, int *num_seqs, FILE *fp)
{
  char buffer[BUFFER_SIZE];
  char *cp = buffer;
  int c;
  int n = 0;
  int new = 1;
  int len = 0;
  int err = GFF_FASTA;
  Seq_list *seq_list;
  Seq_list *seq_list_p;
  Seq *seq = malloc(sizeof(Seq));
  
  //  while ((c = fgetc(fp)) != EOF)
  while (err == GFF_FASTA)
    {
      c = fgetc(fp);
      if (c == '\n')
	{
	  if (new == 1)
	    { // end of ID
	      if (c == EOF) 
		err = GFF_EOF;
	      new = 0;
	      *cp = '\0';
	      printf(":%s",buffer);
	      seq->ID = malloc(sizeof(char)*(len+1));
	      strcpy(seq->ID,buffer);
	      //	      strcpy(seq_list->seq->ID,buffer);
	      //printf("%s\n",seq_list->seq->ID);
	      printf(":%s*\n",seq->ID);
	      cp = buffer;
	      len = 0;
	    }
      	}
      //      else if (c == '>' || c == EOF)
      else if (new == 0 && (c == '>' || c == EOF))
	{ // end of seqstr
	  *cp = '\0';
	  seq->seqstr = malloc(sizeof(char)*(len+1));
	  strcpy(seq->seqstr,buffer);
	  printf("%s len:%d\n",seq->seqstr,len);
	  cp = buffer;
	  len = 0;
	  if (c == EOF)
	    err = GFF_EOF;
	  
	  if (n == 0)
	    {
	      seq_list = (Seq_list*)malloc(sizeof(Seq_list));
	      seq_list_p = seq_list;
	    }
	  else
	    {
	  printf("read fasta2a\n");
	      seq_list_p->next = (Seq_list*)malloc(sizeof(Seq_list));
	      seq_list_p = seq_list_p->next;
	    }
	  seq_list_p->seq = malloc(sizeof(Seq));
	  seq_list_p->seq->ID = seq->ID;
	  seq_list_p->seq->seqstr = seq->seqstr;
	  printf("%s %s\n",seq_list_p->seq->ID,seq_list_p->seq->seqstr);
	  seq_list_p->next = NULL;
	  n++;
	  printf("num_fasta%d\n",n);
	  new = 1;
	}
      else if (new == 1 && (c == '>' || c == EOF))
	{
	  // TODO FREE UP MEMORY
	  //	  err = GFF_BAD_FASTA;
	  printf("read_fasta3a: BAD DATA\n");
	  free(seq->ID);
	  err = GFF_EOF;
	  new = 0;
	}
      else
	{
	  *cp++ = c;
	  len++;
	  printf("%c",c);
	}
    }

  // if we have a GFF_EOF and it isn't at the start it's an error
  if (err == GFF_EOF && new == 0)
    {
    err = GFF_BAD_FASTA;
    seq_list_p = seq_list;
    while (seq_list_p != NULL)
      {
	// free up any malloced memory
	/*
	free(seq_list_p->seq->ID);
	free(seq_list_p->seq->seqstr);
	free(seq_list_p->seq);
	*/
	// TODO possible problem here if seq not fully malloced!
	  printf("read_fasta3b: BAD DATA\n");
	gff_free_seq(seq_list_p->seq);
	  printf("read_fasta3c: BAD DATA\n");
      }
    }
  else // copy linked list to array
    {
      *num_seqs = n;
      //      seqs = malloc(sizeof(Seq*)*n);
      seq_list_p = seq_list;
      //      Seq *seq_p = (Seq*)seqs;
      Seq *seq_p = *seqs;
      len = 0;
      while (seq_list_p != NULL)
	{
	  /*
	  seq_p = seq_list_p->seq;
    printf("%s\n",seq_p->ID);
	  seq_p++;
	  */
    printf("%s\n",seq_list_p->seq->ID);
	  seqs[len++] = seq_list_p->seq;
	  seq_list_p = seq_list_p->next;
	}
    }

  printf("Print Seq ID's %d\n",n);
  for (len=0; len<n; len++)
    printf("%s\n",seqs[len]->ID);
    
  printf("Free up FASTA list\n");
  // free up linked list
  Seq_list *seq_list_p_next;
  seq_list_p = seq_list;
  while (seq_list_p != NULL)
    {
      seq_list_p_next = seq_list_p->next;
      free(seq_list_p);
      seq_list_p = seq_list_p_next;  
    }

  free(seq);
  
  return err;
}


int gff_read_record(Feature *feature, FILE *fp)
{
  int c;
  char buffer[BUFFER_SIZE];
  char *cp = buffer;
  char **columns;
  int x = 0;
  int err = GFF_SUCCESS;
  
  if ((c = fgetc(fp)) == EOF)
    return GFF_EOF;

  printf("gff_read_record:0\n");
  // first check for pragmas and comments
  if (c  == '#')
    {
      //  printf("gff_read_record:1\n");
      if ((c = fgetc(fp)) == '#')
	{
	  // we have a pragma
	  //while ((c = fgetc(fp)) != '\n')
	  //if(c == EOF) return GFF_EOF;
	  err = gff_read_pragma(fp);
	  return err;
	}
      else 
	{
	  // we have a comment
	  while ((c = fgetc(fp)) != '\n')
	    if(c == EOF) return GFF_EOF;
	  return GFF_COMMENT;
	}
    }
  // check to see if the rest of the file is fasta data
  else if (c == '>')
    {
      return GFF_FASTA;
    }
	  
  //printf("gff_read_record:2\n");
  // read each line exiting with GFF_EOF if EOF
  while(c != EOF)
    {
      *cp++ = (char)c;
      if (c  == '\n')
	// then return the Feature
	  {
	    //printf("gff_read_record:3\n");
	    if ((columns = gff_read_columns(buffer)) == NULL)
	      return GFF_BAD_COLUMN;

	    //	    for (x=0; x<9; x++)
	    //	    printf("col %d %s*\n",x,columns[x]);

	    /*
	    // col 1 seqid
	    if (columns[0][0] == '.')
	      {
	      feature->seqid = NULL;
	      free(columns[0]);
	      }
	    else
	      feature->seqid = columns[0];

	    // col2 source
	    if (columns[1][0] == '.')
	      {
	      feature->source = NULL;
	      free(columns[1]);
	      }
	    else
	      feature->source = columns[1];

	    // col3 type
	    if (columns[2][0] == '.')
	      {
	      feature->type = NULL;
	      free(columns[2]);
	      }
	    else
	      feature->type = columns[2];
	    */

	    // col 1 seqid
	    feature->seqid = columns[0];

	    // col2 source
	    feature->source = columns[1];

	    // col3 type
	    feature->type = columns[2];

	    // cols 4,5,7,8 locs(start,end,strand,phase) 
	    feature->locs = malloc(sizeof(SeqLoc*));
	    feature->locs[0] = malloc(sizeof(SeqLoc));
	    int *ip = NULL;

	    if (columns[3][0] == '.')
	      feature->locs[0]->start = NULL;
	    else
	      {
		ip = malloc(sizeof(int));
		sscanf(columns[3],"%d",ip);
		feature->locs[0]->start = ip;
	      }
	    if (columns[4][0] == '.')
	      feature->locs[0]->end = NULL;
	    else
	      {
		ip = malloc(sizeof(int));
		sscanf(columns[4],"%d",ip);
		feature->locs[0]->end = ip;
	      }
	    if (columns[6][0] == '.')
	      feature->locs[0]->strand = NULL;
	    else
	      {
		ip = malloc(sizeof(int));
		feature->locs[0]->strand = ip;
		if (columns[6][0]=='+')
		  *feature->locs[0]->strand = 1; 
		else if (columns[6][0]=='-')
		  *feature->locs[0]->strand = -1; 
		else if (columns[6][0]=='?')
		  *feature->locs[0]->strand = 0;
		else
		  {
		  feature->locs[0]->strand = NULL;
		  free(ip);
		  }
	      }
	    if (columns[7][0] == '.')
	      feature->locs[0]->phase = NULL;
	    else
	      {
	      ip = malloc(sizeof(int));
	      sscanf(columns[7],"%d",ip);
	      feature->locs[0]->phase = ip;
	      }

	    // col 6 score
	    float *f;
	    if (columns[5][0] == '.')
	      feature->score = NULL;
	    else
	      {
		f = malloc(sizeof(float));
		sscanf(columns[5],"%f",f);
		feature->score = f;
	      }
		
	    // col 9 attributes
	    feature->ID = NULL;
	    feature->Name = NULL;
	    feature->parents = NULL;
	    feature->attributes = NULL;
	    char attr_buffer[BUFFER_SIZE];
	    char *colp = columns[8];
	    char attr_name[BUFFER_SIZE];
	    char **attr_vals;
	    int num_vals[BUFFER_SIZE];
	    Attribute *attributes[BUFFER_SIZE];
	    feature->num_parents = 0;
	    feature->num_attributes = 0;
	    //	    while (*colp != '\0')
	    while (*colp != '\0' && *colp != EOF && *colp != '\n')
	      {
		// return an error for bad syntax
		if (*colp == ',' || *colp == ';') return GFF_BAD_COLUMN;

		// for each attribute
		int attr_len = 0;
		char *attrp = attr_buffer;
		//		while (*colp != ';' && *colp != '\0')

		while (*colp != ';' && *colp != '\0' && *colp != '\n'&& *colp != EOF )
		  {
		    *attrp++ = *colp++;
		    attr_len++;
		  }
		attr_buffer[attr_len] = '\0';
		if (strncmp("ID=",attr_buffer,3)==0)
		  {
		    feature->ID = malloc((sizeof(char))*(attr_len-2));
		    strcpy(feature->ID,&attr_buffer[3]);
		  }
		else if (strncmp("Name=",attr_buffer,5)==0)
		  {
		    feature->Name = malloc(sizeof(char)*(attr_len-4));
		    strcpy(feature->Name,&attr_buffer[5]);
		  }
		else if (strncmp("Parent=",attr_buffer,7)==0)
		  {
		    printf("Parent:%s* <%d>\n",colp,colp);
		    feature->parents = gff_split_attribute(&attr_buffer[7],&(feature->num_parents));
		  }
		else
		  {
		    int pos_eq = 0;
		    attrp = attr_buffer;
		    printf("ATTR:%d=%s*\n",feature->num_attributes,attrp);
		    while ((*attrp != '=') && (*attrp != '\0'))
		      {
			attr_name[pos_eq++] = *attrp++;
		      }
		    attr_name[pos_eq] = '\0';
		    attrp++;
		    printf("attr:%d=%s\n",feature->num_attributes,attr_name);
		    attributes[feature->num_attributes] = malloc(sizeof(Attribute));
		    attributes[feature->num_attributes]->tag = malloc(sizeof(char)*(pos_eq+1));
		    //		    attributes[feature->num_attributes]->tag = malloc(sizeof(char)*(pos_eq));
		    strcpy(attributes[feature->num_attributes]->tag,attr_name);
		    //printf("Col8left1:%s* <%d>\n",colp,colp);
		    attributes[feature->num_attributes]->vals = gff_split_attribute(attrp,&attributes[feature->num_attributes]->num_vals);
		    //printf("Col8left2:%s* <%d>\n",colp,colp);
		    //printf("num_vals=%d\n",attributes[feature->num_attributes]->num_vals);
		    feature->num_attributes++;
		  }
		//		colp++;
		if (*colp != '\0') colp++;
		//printf("Col8left3:%s* <%d>\n",colp,colp);
	      }
	    
	    //create attributes
	    feature->attributes = malloc(sizeof(Attribute*)*feature->num_attributes);
	    for (x=0; x<feature->num_attributes; x++)
	      {
		feature->attributes[x] = attributes[x];
		//printf("%s=%s %d\n",feature->attributes[x]->tag,feature->attributes[x]->vals[0],feature->attributes[x]->num_vals);
	      }

	    // free up colums that we don't need
	    for (x=3; x<9; x++)
	      free(columns[x]);

	    return err;
	  }
      c = fgetc(fp);
    }

  printf("gff_read_record:4\n");
  return GFF_EOF;
}

int gff_read_pragma(FILE *fp)
{
  char buffer[BUFFER_SIZE];
  char *cp = buffer;
  int err = GFF_EOF;

  while((*cp = fgetc(fp)) != EOF)
    {
      if (*cp == '\n') 
	{
	  *cp = '\0';
	  printf("PRAGMA=%s*\n",buffer);
	  if (strcmp(buffer,"gff-version 3") == 0)
	    return GFF_DOC;
	  // process PRAGMA
	  return GFF_PRAGMA;
	}
      cp++;
    }

  return err;
}

char **gff_split_attribute(char *attr, int *n)
{
  char *values[BUFFER_SIZE];
  char **ret_values;
  char *attr_start = attr;
  char *attr_end = attr;
  int attr_len = 0;
  int num_values = 0;
  
  printf("gff_split:%s\n",attr);

  while (*attr_end != '\0')
    {
      attr_len++;
      if (*attr_end == ',')
	//      if (*attr_end == ',' || *attr_end == ';')
	{
	  //	  values[num_values] = malloc(sizeof(char)*(attr_len));
	  values[num_values] = malloc(sizeof(char)*(attr_len+1));
	  *attr_end = '\0';
	  strcpy(values[num_values],attr_start);
	  //	  attr_end++;
	  //	  attr_start = attr_end--;
	  attr_start = attr_end+1;
	  attr_len = 0;
	  num_values++;
	}
      attr_end++;
    }
  //  values[num_values] = malloc(sizeof(char)*(attr_len));
  values[num_values] = malloc(sizeof(char)*(attr_len+1));
  *attr_end = '\0';
  strcpy(values[num_values],attr_start);
  num_values++;
  *n = num_values;

  ret_values = malloc(sizeof(char*)*num_values);
  int x=0;
  for (x=0; x<num_values; x++)
    {
    ret_values[x] = values[x];
    //    printf("gff_split:%s\n",ret_values[x]);
    }

  //printf("gff_split:%s*\n",attr_start);
  //printf("gff_split:%s*\n",attr_end);

  return ret_values;
} 

// splits the record into columns, returning NULL on failure
char **gff_read_columns(char *record)
{
  char **columns;
  char *rp = record;
  int x = 0;
  int n = 0;
  char buffer[BUFFER_SIZE];
  char *cp = buffer;

  //printf("gff_read_columns:1\n");
  columns = malloc(sizeof(char*)*9);

  for (x=0; x<9; x++)
    {
      //printf("gff_read_columns:2 ");
      n = 0;
      cp = buffer;
      while (*rp != '\t' && *rp != '\n' && *rp != EOF)
	{
	  *cp++ = *rp++;
	  //printf("%c",*rp); printf("%c",*cp); printf("%c",buffer[n]);
	  n++;
	}
      *cp = '\0';  
      //printf("\ngff_read_columns:3 %s\n",buffer);
      columns[x] = malloc(sizeof(char)*(n+1));
      strcpy(columns[x],buffer);

      // exit clause
      if (*rp == '\n' || *rp == EOF)
	if (x != 8) 
	  {
	    x++; // needed to free last column added
	    break;
	  }
	else 
	  return columns;

      rp++;
    }

  // remember to free up memory if there is a fail
  for (n=0; n<x; n++)
    {
      //      printf("read_column failure freeing column %d\n",n);
      free(columns[n]);
    }
  free(columns);

  return NULL;
}

int gff_free_feature(Feature *fp, int no_locs)
{
  int i,j;
  //printf("gff_free_feature1\n",fp->num_locs);

  //  gff_print_feature(fp);

  free(fp->seqid);
  //printf("gff_free_feature2\n",fp->num_locs);
  free(fp->source);
  //  printf("gff_free_feature3\n",fp->num_locs);
  free(fp->type);
  //  printf("num_locs=%d\n",fp->num_locs);
  if (!no_locs)
    {
      SeqLoc **seqloc  = fp->locs;
      for (i=0; i<fp->num_locs; i++)
	{
	  free((*seqloc)->start);
	  free((*seqloc)->end);
	  free((*seqloc)->strand);
	  free((*seqloc)->phase);
	  free(*seqloc++);
	}
    }
  free(fp->locs);
  free(fp->score);
  free(fp->ID);
  free(fp->Name);
  char **parent = fp->parents;
  for (i=0; i<fp->num_parents; i++)
    {
      //printf("attr=%d\n",i);
      //free(*parent++);
      free(fp->parents[i]);
    }
  free(fp->parents);
  Attribute **attr = fp->attributes;
  for (i=0; i<fp->num_attributes; i++)
    {
      //printf("attr=%d\n",i);
      free((*attr)->tag);
      char **valp = (*attr)->vals;
      for(j=0; j<((*attr)->num_vals); j++)
	{
	  //printf("vals=%d\n",j);
	  //free(*valp++);
	  free((*attr)->vals[j]);
	}
      free((*attr)->vals);
      free(*attr++);
    }
  free(fp->attributes);
  
  //gff_print_feature(*fp);
  free(fp);
}

// free up an individual FASTA sequence
int gff_free_seq(Seq *s)
{
  free(s->ID);
  free(s->seqstr);
  free(s);

  return 1;
}

// frees up all the memory grabbed by the GFFDoc
int gff_free(GFFDoc *gffdoc)
{
  int x = 0;
  Feature **fp = gffdoc->features;
  Seq **sp = gffdoc->seqs;
  int i,j;

  for (x=0; x<gffdoc->num_features; x++)
    {
      // TO DO REPLACE 0 WITH NUM_LOCS
      gff_free_feature(*fp,0);
      fp++;
      printf("Feature %d free\n",x);
    }

  for (x=0; x<gffdoc->num_seqs; x++)
    {
      gff_free_seq(*sp);
      sp++;
      printf("Seq %d free\n",x);
    }

  free(gffdoc->features);
  free(gffdoc->seqs);
  free(gffdoc);

  return GFF_SUCCESS;
}


void gff_print_feature(Feature *f)
{
  int attr = 0;
  int x,loc = 0;

  printf("print feature num_locs=%d\n",f->num_locs);
  for (loc=0; loc<f->num_locs; loc++)
    {
      if (f->seqid == NULL)
	printf(".\t");
      else
	printf("%s\t",f->seqid);

      if (f->source == NULL)
	printf(".\t");
      else
	printf("%s\t",f->source);

      if (f->type == NULL)
	printf(".\t");
      else
	printf("%s\t",f->type);

      if (f->locs[loc]->start == NULL)
	printf(".\t");
      else
	printf("%d\t",*f->locs[loc]->start);

      if (f->locs[loc]->end == NULL)
	printf(".\t");
      else
	printf("%d\t",*f->locs[loc]->end);
      
      if (f->score == NULL)
	printf(".\t");
      else
	printf("%f\t",*f->score);
      
      if (f->locs[loc]->strand == NULL)
	printf(".\t");
      else
	if (*f->locs[loc]->strand == 1)
	  printf("+\t");
	else if (*f->locs[loc]->strand == 0)
	  printf("?\t");
	else if (*f->locs[loc]->strand == -1)
	  printf("-\t");

      if (f->locs[loc]->phase == NULL)
	printf(".\t");
      else
	printf("%d\t",*f->locs[loc]->phase);

      if (f->ID != NULL)
	{
	  //      printf("ID%d=%s*",attr,f->ID);
	  printf("ID=%s*",f->ID);
	  attr++;
	}
      if (f->Name != NULL)
	{
	  if (attr > 0) printf(";"); 
	  printf("Name=%s",f->Name);
	  attr++;
	}
      for (x=0; x<f->num_parents; x++)
	{
	  if (attr > 0 && x == 0)
	    printf(";Parent=");
	  else if (x == 0) 
	    printf("Parent=");
	  if (x > 0) printf(",");
	  printf("%s*",f->parents[x]);
	  attr++;
	}
      for (x=0; x<f->num_attributes; x++)
	{
	  if (attr > 0) printf(";");
	  //      printf("Attrs%d:",f->num_attributes);
	  printf("%s=",f->attributes[x]->tag);
	  int i = 0;
	  for (i=0; i<f->attributes[x]->num_vals; i++)
	    {
	      if (i > 0)
		printf(",");
	      printf("%s*",f->attributes[x]->vals[i]);
	    }
	  attr++;
	}
      printf("\n");
      
    }
  return;
}

void gff_print(GFFDoc *gffdoc)
{
  int x = 0;
  Feature **fp = gffdoc->features;

  printf("Print file\n");
  for (x=0; x<gffdoc->num_features; x++)
    {
      gff_print_feature(*fp++);
    } 
  
  // print fasta
  printf("num_seqs=%d\n",gffdoc->num_seqs);
  Seq **seq_p = gffdoc->seqs;
  for (x=0; x<gffdoc->num_seqs; x++)
    {
      //      printf("%s : %s\n",gffdoc->seqs[x]->ID,gffdoc->seqs[x]->seqstr);
      printf("%s : %s\n",(*seq_p)->ID,(*seq_p)->seqstr);
      seq_p++;
    } 

  return;
}
