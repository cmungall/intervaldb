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
  GFFDoc *gff = malloc(sizeof(GFFDoc));
  
  err = gff_get_doc(gff,"t/gfftab.dat");

  if (err) // positive value so ok
    {
      printf("Print file\n");
      for (x=0; x<gff->num_features; x++)
	gff_print_feature(gff->features[x]);
    }
  err = gff_free(gff);
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

  // allocate enough memory for more records than we need
  gffdoc->features = malloc(sizeof(Feature*)*n);
  err = gff_read_file(gffdoc->features, &gffdoc->num_features, file);

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

int gff_read_file(Feature **features,int *n,char *file)
{
  FILE *fp;
  int x = 0;
  Feature *feature;
  int err = GFF_SUCCESS;

  // open file and count lines
  if (!(fp = fopen(file,"r"))) return GFF_NO_FILE;
  
  // first check that we are dealing with a GFF 3 file
  feature = malloc(sizeof(Feature)); 
  if ((err = gff_read_record(feature,fp)) != GFF_DOC)
    {
      free(feature);
      return GFF_FAIL;
    }

  while (err != GFF_EOF)
    {
      printf("gff_read_file:1 x=%d\n",x);
      feature = malloc(sizeof(Feature)); 
      err = gff_read_record(feature,fp);
      switch(err) {
      case GFF_SUCCESS: 
	printf("line %d\n",x);
	features[x++] = feature;
	break;
      case GFF_PRAGMA: 
	printf("line %d is a PRAGMA\n",x);
	free(feature);
	break;
      case GFF_COMMENT: 
	printf("line %d is a Comment\n",x);
      }
    }

  free(feature);
  *n = x;

  printf("gff_read_file:2\n");
  fclose(fp);
  printf("gff_read_file:3\n");

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
      printf("gff_read_record:1\n");
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
	  
  printf("gff_read_record:2\n");
  // read each line exiting with GFF_EOF if EOF
  while(c != EOF)
    {
      *cp++ = (char)c;
      if (c  == '\n')
	// then return the Feature
	  {
	    printf("gff_read_record:3\n");
	    if (columns = gff_read_columns(buffer))
	      {
		for (x=0; x<9; x++)
		  printf("col %d %s*\n",x,columns[x]);
		// col 1 seqid
		feature->seqid = columns[0];

		// col2 source
		feature->source = columns[1];

		// col3 type
		feature->type = columns[2];

		// cols 4,5,7,8 locs(start,end,strand,phase) 
		feature->locs = malloc(sizeof(SeqLoc*));
		feature->locs[0] = malloc(sizeof(SeqLoc));
		int i;
		sscanf(columns[3],"%d",&i);
		feature->locs[0]->start = i;
		sscanf(columns[4],"%d",&i);
		feature->locs[0]->end = i;
		if (columns[6][0]=='+')
		  feature->locs[0]->strand = 1; 
		else if (columns[6][0]=='-')
		  feature->locs[0]->strand = -1; 
		else if (columns[6][0]=='?')
		  feature->locs[0]->strand = 0; 
		else 
		  feature->locs[0]->strand = 0; 
		if(sscanf(columns[7],"%d",&i))
		  feature->locs[0]->phase = i;
		else
		  feature->locs[0]->phase = -1;

		// col 6 score
		int f;
		if (sscanf(columns[5],"%f",&f))
		  feature->score = f;
		else
		  feature->score = 0;

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
		while (*colp != '\0')
		  {
		    // for each attribute
		    int attr_len = 0;
		    char *attrp = attr_buffer;
		    while (*colp != ';' && *colp != '\0')
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
			feature->parents = gff_split_attribute(&attr_buffer[7],&(feature->num_parents));
			//err = gff_split_attribute(feature->parents, &attr_buffer[7],&(feature->num_parents));
		      }
		    else
		      {
			int pos_eq = 0;
		        attrp = attr_buffer;
			while (*attrp != '=' && *attrp != '\0')
			  {
			    attr_name[pos_eq] = *attrp++;
			    pos_eq++;
			  }
			attr_name[pos_eq] = '\0';
			attrp++;
			printf("attr:%s\n",attr_name);
			attributes[feature->num_attributes] = malloc(sizeof(Attribute));
			attributes[feature->num_attributes]->tag = malloc(sizeof(char)*(pos_eq));
			strcpy(attributes[feature->num_attributes]->tag,attr_name);
			attributes[feature->num_attributes]->vals = gff_split_attribute(attrp,&attributes[feature->num_attributes]->num_vals);
			// ERROR IS HERE
			//			err = gff_split_attribute(attributes[feature->num_attributes]->vals, attrp, &attributes[feature->num_attributes]->num_vals);
			printf("num_vals=%d\n",attributes[feature->num_attributes]->num_vals);
			feature->num_attributes++;
		      }
		    colp++;
		  }

		//create attributes
		feature->attributes = malloc(sizeof(Attribute*)*feature->num_attributes);
		for (x=0; x<feature->num_attributes; x++)
		  {
		    feature->attributes[x] = attributes[x];
		    printf("%s=%s %d\n",feature->attributes[x]->tag,feature->attributes[x]->vals[0],feature->attributes[x]->num_vals);
		  }

		// free up colums that we don't need
		free(columns[3]);
		free(columns[4]);
		free(columns[5]);
		free(columns[7]);
		free(columns[8]);

		return err;
	      }
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
  
  while (*attr_end != '\0')
    {
      attr_len++;
      if (*attr_end == ',')
	{
	  values[num_values] = malloc(sizeof(char)*(attr_len));
	  *attr_end = '\0';
	  strcpy(values[num_values],attr_start);
	  attr_end++;
	  attr_start = attr_end--;
	  attr_len = 0;
	  num_values++;
	}
      attr_end++;
    }
  values[num_values] = malloc(sizeof(char)*(attr_len));
  *attr_end = '\0';
  strcpy(values[num_values],attr_start);
  num_values++;
  *n = num_values;

  ret_values = malloc(sizeof(char*)*num_values);
  int x=0;
  for (x=0; x<num_values; x++)
    {
    ret_values[x] = values[x];
    printf("%s\n",ret_values[x]);
    }
  return ret_values;
} 

char **gff_read_columns(char *record)
{
  char **columns;
  char *rp = record;
  int x = 0;
  int n = 0;
  char buffer[BUFFER_SIZE];
  char *cp = buffer;

  printf("gff_read_columns:1\n");
  columns = malloc(sizeof(char*)*9);

  for (x=0; x<9; x++)
    {
      printf("gff_read_columns:2 ");
      n = 0;
      cp = buffer;
      while (*rp != '\t' && *rp != '\n')
	{
	  //	  *cp = (char)*rp;
	  buffer[n] = *rp;
	  //       	  printf("%c",*rp); printf("%c",*cp); printf("%c",buffer[n]);
	  n++; cp++; rp++;
	}
      *cp = '\0';  
      printf("\ngff_read_columns:3 %s\n",buffer);
      columns[x] = malloc(sizeof(char)*(n+1));
      strcpy(columns[x],buffer);
      if (*rp++ == '\n')
	return columns;
    }

  return columns;
}

int gff_free(GFFDoc *gffdoc)
{
  return GFF_SUCCESS;
}


void gff_print_feature(Feature *f)
{
  int attr = 0;
  printf("%s\t%s\t%s\t%d\t%d\t%f\t%d\t%d\t",f->seqid,f->source,f->type,f->locs[0]->start,f->locs[0]->end,f->score,f->locs[0]->strand,f->locs[0]->phase);
  if (f->ID != NULL)
    {
      printf("ID=%s",f->ID);
      attr++;
    }
  if (f->Name != NULL)
    {
      if (attr > 0) printf(";"); 
      printf("Name=%s",f->Name);
      attr++;
    }
  int x=0;
  for (x=0; x<f->num_parents; x++)
    {
      if (attr > 0 && x == 0)
	printf(";Parent=");
      else if (x == 0) 
	printf("Parent=");
      if (x > 0) printf(",");
      printf("%s",f->parents[x]);
      attr++;
    }
  for (x=0; x<f->num_attributes; x++)
    {
      if (attr > 0) printf(";");
      printf("%s=",f->attributes[x]->tag);
      int i = 0;
      for (i=0; i<f->attributes[x]->num_vals; i++)
	{
	  if (i > 0)
	    printf(",");
	  printf("%s",f->attributes[x]->vals[i]);
	}
      attr++;
    }
  printf("\n");
  return;
}
