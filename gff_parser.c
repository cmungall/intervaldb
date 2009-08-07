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
//#define GFF_EOF -1

main()
{
  GFFDoc *gff;
  //  gff = gff_get_doc("gff3.dat");
  gff = gff_get_doc("gfftab.dat");
}

GFFDoc *gff_get_doc(char *file)
{
  GFFDoc *gffdoc = NULL;
  int x = 0;

  // create the GFFDoc failing with 0
  if (gffdoc = malloc(sizeof(GFFDoc)))
    {
      // count the number of lines and create an array of Features
      int n = gff_num_lines(file);
      printf("number of lines = %d\n",n);
  
      int r = gff_read_lines(n,file);
      printf("number of records = %d\n",r);

      gffdoc->features = gff_read_file(r,file);

      printf("Print file\n");
      for (x=0; x<r; x++)
      	gff_print_feature(gffdoc->features[x]);
    }

  return gffdoc;
}

int gff_num_lines(char* file)
{
  FILE *fp;
  int n = 0;
  int c;

  // open file and count lines
  if (!(fp = fopen(file,"r"))) printf("File not opened\n");

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
  if (!(fp = fopen(file,"r"))) printf("File not opened\n");

  while (i<n && fscanf(fp,"%s\n",&line))
    {
      i++;
      // check line is not a comment or empty
      if (line[0] != '#' && line[0] != '\n') 
	{
	r++;
	printf("%s\n",line);
	}
    }

  fclose(fp);
  return r;
}

Feature **gff_read_file(int r,char *file)
{
  FILE *fp;
  int x = 0;
  Feature **features;
  Feature *feature;

  // create an array of Features
  if (features = malloc((sizeof(Feature)*r)))
    {
      printf("gff_read_file:0 r=%d\n",r);
      // open file and count lines
      if (!(fp = fopen(file,"r"))) printf("File not opened\n");

      //  for (x=0; x<r; x++)
      while (x<r)
      {
	printf("gff_read_file:1 x=%d\n",x);
	if ((feature = gff_read_record(fp)) != NULL)
	  { 
	    printf("line %d\n",x);
	    features[x++] = feature;
	  }
	else
	  return NULL;
      }
    }

  printf("gff_read_file:2\n");
  fclose(fp);
  printf("gff_read_file:3\n");
  return features;
}

Feature *gff_read_record(FILE *fp)
{
  Feature *feature = NULL;
  int c;
  char buffer[BUFFER_SIZE];
  char *cp = buffer;
  char **columns;
  int x = 0;
  
  if ((c = fgetc(fp)) == EOF)
    //    return GFF_EOF;
    return NULL;

  printf("gff_read_record:0\n");
  // first check for pragmas and comments
  if (c  == '#')
    {
      printf("gff_read_record:1\n");
      if ((c = fgetc(fp)) == '#')
	{
	  // we have a pragma
	  while ((c = fgetc(fp)) != '\n')
	    if(c == EOF) return NULL;
	  feature = gff_read_record(fp);
	  return feature;
	}
      else 
	{
	  // we have a comment
	  while ((c = fgetc(fp)) != '\n')
	    if(c == EOF) return NULL;
	  feature = gff_read_record(fp);
	  return feature;
	}
    }
	  
  printf("gff_read_record:2\n");
  // read each line exiting NULL if EOF
  while(c != EOF)
    {
      *cp++ = (char)c;
      if (c  == '\n')
	// then return the Feature
	if(feature = (Feature*)malloc(sizeof(Feature)))
	  {
	    printf("gff_read_record:3\n");
	    if (columns = gff_read_columns(buffer))
	      {
		for (x=0; x<9; x++)
		  printf("col %d %s*\n",x,columns[x]);
		feature->seqid = columns[0];
		feature->source = columns[1];
		feature->type = columns[2];
		feature->locs = malloc(sizeof(SeqLoc*));
		feature->locs[0] = malloc(sizeof(SeqLoc));
		int i;
		sscanf("1000","%d",&i);
		sscanf(columns[3],"%d",&i);
		feature->locs[0]->start = i;
		sscanf(columns[4],"%d",&i);
		feature->locs[0]->end = i;
		int f;
		if (sscanf(columns[5],"%f",&f))
		  feature->score = f;
		else
		  feature->score = 0;
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

		// now attributes
		feature->ID = NULL;
		feature->Name = NULL;
		feature->parents = NULL;
		feature->attributes = NULL;
		char attr_buffer[BUFFER_SIZE];
		char *colp = columns[8];
		//		char *parents[BUFFER_SIZE];
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
			//parents[feature->num_parents] = malloc(sizeof(char)*(attr_len-6));
			//strcpy(parents[feature->num_parents],&attr_buffer[7]);
			//feature->num_parents++;
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
			printf("num_vals=%d\n",attributes[feature->num_attributes]->num_vals);
			feature->num_attributes++;
		      }
		    colp++;
		  }

		//create parents and attributes
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

		return feature;
	      }
	  }
      c = fgetc(fp);
    }

  printf("gff_read_record:4\n");
  return NULL;
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

void gff_print_feature(Feature *f)
{
  int attr = 0;
  printf("%s\t%s\t%s\t%d\t%d\t%f\t%d\t%d\t",f->seqid,f->source,f->type,f->locs[0]->start,f->locs[0]->end,f->score,f->locs[0]->strand,f->locs[0]->phase);
  //printf("%s\t%s\t%s\n",f->seqid,f->source,f->type);
  //  printf("%s\n",f->seqid);
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
