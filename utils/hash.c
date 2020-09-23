/* ------------------------------------------------------------------ *
 * Permission to copy and/or modify all or part of this work is       *
 *   granted, provided that this NO WARRANTY and this copyright	      *
 *   notice are retained verbatim and are displayed conspicuously. If *
 *   anyone needs other permissions that aren't covered by the above, *
 *   please contact the author					      *
 *      							      *
 * NO WARRANTY: THIS WORK IS PROVIDED ON AN "AS IS" BASIS.  THE       *
 *   AUTHOR PROVIDES NO WARRANTY WHATSOEVER, EITHER EXPRESS OR        *
 *   IMPLIED, REGARDING THE WORK, INCLUDING WARRANTIES WITH THE WORK, *
 *   INCLUDING WARRANTIES WITH RESPECT TO ITS MERCHANTABILITY OR      *
 *   FITNESS FOR ANY PARTICULAR PURPOSE.			      *
 *								      *
 * Philip Lijnzaad, lijnzaad@embl-heidelberg.de			      *
 * ------------------------------------------------------------------ */

/* module for doing simple hash table stuff. Both keys and values are */
/* either fixed length, or null terminated strings; this is specified */
/* when setting up the table, with Hcreate. A keysize of zero implies a */
/* 0-terminated string; keysize == -1 implies a pointer (or comparable) */ 

/* this module should be re-written to use linear open adressing, instead
   of linked lists dangling from buckets. Probably faster, more
   economical, esp. if quadratic open adressing is used: if bucket B[i] is
   full, search (B[i + j*j])%n and B[i - j*j]%n for 1 <= j <= (n-1)/2. 
   N, the nr of buckets, must be a prime of form 4*k - 3. The advantage of
   this scheme is that 'clustering' of occupied buckets is avoided better.
   You still can move the currently found one to the 'right' bucket to
   improve the hit chance, by exchanging. */

#include <string.h> 
#include <limits.h>			/* for definition of USHRT_MAX */
#include <ctype.h> 
#include <stdlib.h>

#include "utils.h"
#include "alloc.h" 
#include "hash.h"

static struct  {
  int prime_table_size;
  ushort * prime_table;			/* ushort because may get big */
  int totsize, leftover;		/* for allocation of prime_table */
  int cmp_key_size;			/* has to be static ! */
} local = { 1, NULL, 0,0, -1 };

static int ptr_cmp(void *a, void *b) {
  return ( (2*(a>b)-1)*(a!=b)  );
}

typedef int (*function_pointer)();	/* generic function pointer */

/* static function_pointer diff_func[3]= { memcmp, strcmp, ptr_cmp }; */


#define DIFF_FUNC(KS) ((KS)==0?(&strcmp):((KS)<0?(&ptr_cmp):(&memcmp)))

/* local prototypes: */
static int prime_above(int n);
static int lookup_prime(int n);
static void write_item(FILE *file, const void * item, char * fmt, int size);
static void skip_comment(FILE *file);
static const void * read_item(FILE *file, char * fmt, int size);
int fgetz(char * string, FILE* file);

#if ! defined __GNUC__
#define inline				/* mask gcc's inline directive */
#endif 

#define SEPARATOR '\0'			/* as usual */
#define MAX_ITEM_SIZE 512		/* max size of a key or value */
#define BLKSIZ 512 			/* by which to step up allocation */
					/* of primes */

static inline unsigned int hash(uint tablesize, 
				const char * data, int len) {
  /* adapted from PJ Weinbergers hash function for 0-terminated strings, */
  /* as in Aho, Sethi, Ulmann, p. 434-437, as found in gawk/array.c
   */
  unsigned long h = 0, g;
  const char * s;
  int n;

#define LOOP(TEST) while (TEST) { 		\
		     h = (h << 4) + *s++;	\
		     g = (h & 0xf0000000);	\
		     if (g) {			\
		       h = h ^ (g >> 24);	\
		       h = h ^ g;		\
		     }				\
		   }

  if (len<0) {				/* direct pointer: pretend we */
    s= (const char*)&data;		/* have a pointer to it */
    n=sizeof(void*);
  }
  else {
    s= (char*)data;
    n=len;
  }

  if(len>0)				/* fixed length record */
    { LOOP(n--)	}
  else					/* 0-terminated string */
    { LOOP(*s)	}

  h %= tablesize;
  return h;
} /* hash */

#ifdef FASTHASH
static inline unsigned int hash2(uint tablesize, void * data, int len) {
  /* 2nd hash function, to speed up discarding candidates in bucket when */
  /* list is long (so we don't have to call memcmp on keys time and again) */
  /* This function can be kept simple, since most of the work is done */
  /* already by hash(). Idea nicked from Larry Wall's source code to Perl */

  if (len < 0)
    die("hash2 not implemented for key_size < 0\n");

  char *s = data;
  int n= len;
  int h=0;

  if (len)				/* fixed length record */
    while (n--)
      h = h*5 + *s++;
  else					/* 0-terminated string */
    while (*s)
      h = h*5 + *s++;
    
/* no reducing to tablesize here !! */
  return h;
} /* hash2 */
#endif 


htable_p Hcreate(htable_p init, int wantedsize, int key_size 
		 /* ,  int value_size */ ) {
  /* also serves as Hcopy */
  htable_p ht;
  hnode_p node;
  int size;

  ht=(htable_p)CALLOC(1,sizeof(htable_t));

  if (wantedsize > (int)USHRT_MAX)
    die("Hcreate: requesting table larger than USHRT_MAX\n");

  ht->size=size=prime_above(wantedsize);

  ht->key_size=key_size;
/*  ht->value_size=value_size; */
  ht->difference_test=DIFF_FUNC(key_size); /* strcmp or memcmp; */
  /* you can ignore "pointer type mismatch" message; they're function
     pointers which are impossible to cast appropriately on most compilers */
  ht->table=(hnode_p*)CALLOC(size, sizeof(hnode_p));
  ht->unused=size;

  if(init) {
    if(init->key_size != key_size /* || init->value_size != value_size) */ )
      die("trying to initialize with incompatible hashtable\n");
    
    for(node=Hkeys(init); node; node=node->nexthnode) {
      Hset(ht, node->key)->value=node->value;

      /* for a 'deep copy':
       * 	*Hset(ht, node->key)=memsav(node->value, value_size); 
       */ 
    }
  }
  
  return(ht);
} /* Hcreate */


#ifdef GC
#define STANDARD_FUNC GC_free
#else
#define STANDARD_FUNC _free
#endif /* GC */

void _Hfree(htable_p ht, void(*free_func)(), 
	    const char *file, int line) {
  /* deletes complete htable, using free_func (unless NULL) on list->ptr */
  hnode_p keys, node, next;
  int size=ht->key_size;

  keys=Hkeys(ht);

#define FREEUP(X)\
  if(free_func==(void(*)(void*))STANDARD_FUNC) { /* (jaysus ...) */	 \
    for(node=keys; node; node=node->next)				 \
      if (node->X)			/* may be NULL */		 \
	(*STANDARD_FUNC)(node->X, file, line); /* gets dmw_info right */ \
  }  else				/* custom free function */	 \
    for(node=keys; node; node=node->next)				 \
      (*free_func)(node->X);		/* should cope itself with NULL */

  if(free_func) {			/* free individual data, if present */
    if (size < 0) { 
      FREEUP(key);
    } else {
      FREEUP(value);
    }
  }

  for (node=keys; node; node=next) {
    next=node->next;
    if (size>=0)
      FREE(node->key);
    FREE(node);
  }

  FREE(ht->table);			/* pointer table */
  FREE(ht);				/* the struct itself */
  return;
} /* Hfree */

hnode_p Hkeys(htable_p ht) {
  /* return keys of ht as linked list.  */
  int i, size=ht->size, maxlen=0, len;
  hnode_p *table=ht->table, *prevnext=NULL, keys=NULL, bucket, node;

  if (ht->keys)				/* still valid; may even be sorted */
    return(ht->keys);

  prevnext= &keys;

  for (i=0; i<size; i++) {
    bucket=table[i];
    if ( bucket == NULL )
      continue;
    *prevnext = bucket;			/* modify previous end of list */

    for(len=1, node=bucket; node->nexthnode; node=node->nexthnode, len++) 
      node->next=node->nexthnode;
    node->next=NULL;			/* may point to invalid node !!! */
    prevnext = &(node->next);		/* address of new end of the list */

    if(len>maxlen)			/* update max */
      maxlen=len;
  } /* for */

  ht->keys=keys;			/* cache it */
  ht->sorted_keys=NULL;
  ht->maxlen=maxlen;
  return(keys);				/* may be NULL */
} /* Hkeys */

static int Hmemcmp(hnode_p*a, hnode_p*b) { /* fixed length records */
  return memcmp((char*)(*a)->key, (char*)(*b)->key, local.cmp_key_size);
}
static int Hstrcmp(hnode_p*a, hnode_p*b) { /* 0-terminated strings */
  return strcmp((char*)(*a)->key, (char*)(*b)->key);
}
static int Hptr_cmp(hnode_p*a, hnode_p*b) { /* direct pointers */
  return ptr_cmp((*a)->key, (*b)->key);
}

#define COMPARE_FUNC(KS) ((KS)==0?(Hstrcmp):((KS)<0?(Hptr_cmp):(Hmemcmp)))

hnode_p Hsorted_keys(htable_p ht, int (*compare)(void*, void*)) {
  int i, len;
  hnode_p *array, keys, node;
  int (*cmp)();

  if(ht->sorted_keys)
    return ht->sorted_keys;

  keys=Hkeys(ht);

  len=ht->nkeys;

#ifndef NO_ALLOCA
  array=(hnode_p*)alloca(len*sizeof(hnode_p));
#else
  array=(hnode_p*)MALLOC(len*sizeof(hnode_p));
#endif

  i=0;
  for(node=keys; node; node=node->next)
    array[i++]=node;

  if(compare)
    cmp=compare;			/* used that */
  else { 
    cmp=COMPARE_FUNC(ht->key_size);	/* choose based on key_length */
    local.cmp_key_size=ht->key_size;
  }
  qsort(array, len, sizeof(hnode_p), cmp);

  for (i=0; i<=len-2; i++) 
    array[i]->next=array[i+1];
  array[len-1]->next=NULL;
  keys=array[0];
  ht->keys=keys;
  ht->sorted_keys=keys;

#ifdef NO_ALLOCA
  FREE(array);
#endif 

  return keys;
} /* Hsorted_keys */

#ifdef FASTHASH				/* secondary index for speed */
#  define FASTHASH_EXTRA0    , h2
#  define FASTHASH_EXTRA1   h2 = hash2(ht->size, key, ks)
#  define FASTHASH_EXTRA2   if(node->second_hash != h2) continue
#else
#  define FASTHASH_EXTRA0
#  define FASTHASH_EXTRA1
#  define FASTHASH_EXTRA2
#endif

#define FINDNODE \
  table=ht->table;							  \
  difference=ht->difference_test;					  \
  ks= ht->key_size;							  \
  h = hash(ht->size, key, ks);						  \
  FASTHASH_EXTRA1;							  \
									  \
  prev=NULL;								  \
  for(node=table[h]; node;  node=node->nexthnode) { 			  \
    FASTHASH_EXTRA2;							  \
    if( node->key==key  || (ks>=0 && (*difference)(node->key, key, ks)==0)) \
      break;								  \
    prev=node;								  \
  }

void * _Hdelete(htable_p ht, const void *key, 
		const char * file, int line) {
  /* deletes entry key from htable, returning the value ptr; may need */
  /* freeing itself. Returns NULL if key not found */
  int ks;
  uint h FASTHASH_EXTRA0 ;

  hnode_p node, prev, *table;
  function_pointer difference;
  void * value;

  FINDNODE;
  if(!node)				/* not found */
    return NULL;
  /* else */
  ht->keys=NULL;			/* not valid anymore, so make NULL */
  ht->sorted_keys=NULL;			/* not valid anymore, so make NULL */
  ht->nkeys--;

  if (prev)				/* not first of list */
    prev->nexthnode=node->nexthnode;	/* simple */
  else 
    table[h]=node->nexthnode;	/* first of list */

  ht->unused += (table[h]==NULL);	/* increment if slot now NULL */

  value=node->value;
  _free(node->key, file, line);
  _free(node,  file, line);
  return value;
} /* Hdelete */

hnode_p Hget(htable_p ht, const void *key) {
  int ks;
  uint h FASTHASH_EXTRA0 ;

  hnode_p node, prev, *table;
  function_pointer difference;

  FINDNODE;
  if(!node)				/* not found */
    return NULL;

  if (prev) { 				/* not first of list: make it first */
    ht->keys=NULL;			/* now invalid */
    ht->sorted_keys=NULL;		/* not valid anymore, so make NULL */
    prev->nexthnode=node->nexthnode;
    node->nexthnode=table[h];
    table[h]=node;
  }
  /* else: already first of list */
  return node;
} /* Hget */

hnode_p Hfeel(htable_p ht, const void *key) {
  int ks;
  uint h FASTHASH_EXTRA0 ;
  hnode_p node, *table, prev;
  function_pointer difference;

  FINDNODE;
  if(!node)				/* not found */
    return NULL;
  /* don't reorder: return now */

  return node;
} /* Hfeel */

hnode_p _Hset(htable_p ht, const void *key, 
	       const char * file, int line ) {
  int ks;
  uint h FASTHASH_EXTRA0 ;
  hnode_p node, prev, *table;
  function_pointer difference;

  FINDNODE;
  if(node)	{			/* found: don't have to make one */
    if (prev) {				/* not first of list: make it first */
      ht->keys=NULL;			/* now invalid */
      ht->sorted_keys=NULL;		/* not valid anymore, so make NULL */
      prev->nexthnode=node->nexthnode;
      node->nexthnode=table[h];
      table[h]=node;
    }
    /* else: already first of list */
    return node;
  }
  /* else: not found. Make one, and insert as new head */
  node=_calloc(1, sizeof(hnode_t), file, line);

#ifdef FASTHASH
  node->second_hash = h2;
#endif 
  if (ht->key_size < 0)
    node->key = key;			/* just copy pointer (discarding
					   const qualifier is OK) */
  else
    node->key = _memsav(key, ht->key_size, file, line);

  ht->keys=NULL;
  ht->sorted_keys=NULL;			/* not valid anymore, so make NULL */ 
  ht->nkeys++;
  node->nexthnode=table[h];
  ht->unused -= (table[h]==NULL);	/* decrement it if slot was NULL */
  table[h]=node;

  return node;
} /* Hset */

#define MAYBE_EXTEND_ALLOC(PTR) if(!--local.leftover) {\
PTR=(ushort*)REALLOC(PTR,(local.totsize += BLKSIZ)*sizeof(ushort));local.leftover=BLKSIZ;}

static int prime_above(int n) {
  /* returns first prime larger than n */
  /* calculates primes by Eratosthen's seave */

  int i, j, k, last;

  if(!local.prime_table) {
    local.prime_table=(ushort*)MALLOC(BLKSIZ*sizeof(ushort));
    local.prime_table[0]=3;			/* first prime */
    local.leftover=local.totsize=BLKSIZ;
  }

/*   if (n > USHRT_MAX)			/# already testet in Hcreate! #/
 *     die("requesting prime larger than USHRT_MAX\n");
 */
  if ( (last=local.prime_table[ local.prime_table_size -1 ]) > n) /* already computed */
    return lookup_prime(n) ;

  /* else: compute and extend local.prime_table */

  for (i= last+2; 1; i+=2) {		/* i is the prime candidate checked */
    for (j=0; j < local.prime_table_size; j++) { /* walk local.prime_table */
      k=local.prime_table[j];			/* a possible factor of i */
      if (i % k == 0)			/* not a prime: 2-fold break */
	goto nexti ;			/* try next prime candidate */
      if (SQR(k) > i ) {		/* don't have to look further */
	MAYBE_EXTEND_ALLOC(local.prime_table); /* may need extending */
	local.prime_table[ local.prime_table_size++ ]=i;
	if(i >= n)
	  return i;
	goto nexti;			/* next prime candidate */
      }
    }
  nexti:
    ;
  }
} /* prime_above */

static int lookup_prime(int n) { 
  /* adapted from Numerical Recipes 1st ed., procedure 'locate' */
  int l, u, m;

  /* search by bisection */
  l=0; u=local.prime_table_size;
  while (u-l>1) {
    m=(u+l)>>1;				/* midpoint */
    if (n > local.prime_table[m])
      l=m;
    else
      u=m;
  }
  return local.prime_table[ l +1 ];
} /* lookup_prime */

void Hwrite(FILE*file, htable_p table, char *keyfmt, char *valfmt,
	    int valsize) {
  hnode_p node;
  int keysize=table->key_size;

  if(table->key_size < 0)
    die("Hwrite not implemented for key_size < 0\n");
  for(node=Hkeys(table); node; node=node->next) {
    write_item(file, node->key, keyfmt, keysize);
    write_item(file, node->value, valfmt, valsize);
  }
  return;
} /* Hwrite */

static void write_item(FILE *file, const void * item, char * fmt, int size) {
  int len, i;
  char string[MAX_ITEM_SIZE];

  len = size ? size : strlen( (char*)item);

  if (fmt==NULL || fmt[0]=='\0') {	/* write binary pattern */
    fwrite(item, len, 1, file);
    if (size==0)			/* insert separator character */
      putc(SEPARATOR, file);
    return;
  }

  if (size==0)
    sprintf(string, fmt, (char*)item);	/* item is a char* per definition */
  else {				/* have to find out what it is  */
    switch(size) {			/* pass args to sprintf as if they */
					/* were of prototypical size */
    case 1:
      sprintf(string, fmt, *(char*)item); break;
    case 2:
      sprintf(string, fmt, *(short*)item); break;
    case 4:
      sprintf(string, fmt, *(int*)item); break;
    case 8:
      sprintf(string, fmt, *(double*)item); break;
    default: 
      DIE("item '%0xlXs' has unknown size %d", *(unsigned long*)item, size);
    }
  }
  i=strlen(string);
  if (i > MAX_ITEM_SIZE)
    DIE("length of string after sprintf conversion too large");

  fwrite(string, i, 1, file);
  return;
} /* write_item */

int Hread(FILE*file, htable_p table, char *key_fmt, char *value_fmt,
	   int value_size, 
	  void*(*parse_func)(const char *string) ) {
  int key_size=table->key_size, n, ascii;
  const void *key, *value, *v;
  hnode_p node=NULL, prev_node=NULL;

  assert(table);			/* should have been init-ed */

  if(key_size < 0)
    die("Hread not implemented for key_size < 0\n");

  ascii= key_fmt&&key_fmt[0];		/* i.e. not binary */

  if (key_size > MAX_ITEM_SIZE || value_size > MAX_ITEM_SIZE) 
    DIE("size of key or value too great");

#define CHECK(S) if(!S) {						\
  warn("could not read %s as '%s'\n", #S, S##_fmt);			\
  if (prev_node) {							\
    fprintf(stderr,"(previously read key-value pair was |");		\
    write_item(stderr, prev_node->key, key_fmt, key_size);		\
    fprintf(stderr,"|");						\
    write_item(stderr, prev_node->value, value_fmt, value_size);	\
    fprintf(stderr,"|)\n");						\
  } else warn("(first key-value pair read)\n");				\
  die(""); }

  n=table->nkeys;
  while(1) {
    prev_node=node;
    if (ascii)
      skip_comment(file);
    if (feof(file))
      break;
    key=read_item(file, key_fmt, key_size);
    if (feof(file))
      break;
    CHECK(key);
    if ((node=Hget(table, key))) {
      fprintf(stderr, "key '"); 
      write_item(stderr, node->key, key_fmt, key_size);
      die("' already exists\n");
    }
    node=Hset(table, key);		/* create entry */
    value=read_item(file, value_fmt, value_size);
    CHECK(value);
    if (parse_func != NULL) {
      v=value;
      value= (*parse_func)(v);
      CHECK(value);		 	/* again; maybe error during parse */
      node->value=value;		/* discarding const qualifier OK */
    }
    else
      node->value=memsav(value, value_size);
  } /* while 1 */
  return table->nkeys - n;		/* nr of new entries */
} /* Hread */

static const void * read_item(FILE *file, char * fmt, int size) {
  int i;
  static char string[MAX_ITEM_SIZE];

  if (fmt==NULL || fmt[0]=='\0') {	/* read binary pattern */
    if (size==0) {			/* read until SEPARATOR */
      i=fgetz(string,file);		/* read SEP.-terminated string */
      if (i > MAX_ITEM_SIZE)
	die("length of string read too large\n");
      i=(i>0);				/* convert to boolean */
    }
    else
      i=fread(string, size, 1, file);	/* read fixed nr of bytes */
    if (i!=1)
      return NULL;
    return string;
  }

  /* else */ 				/* read formatted */
  if(fscanf(file, fmt, string)!=1)	/* exactly one item */
    return NULL;

  if (size==0 && strlen(string) > MAX_ITEM_SIZE)
    die("length of string read too large\n");

  /* else: just trust it's an integer or a double or so, can't get too long */
  return string;
} /* read_item */

int fgetz(char * string, FILE* file) {
  /* meant to replace fscanf(file, "%s", string) such that not white */
  /* space, but SEPARATOR makes it stop. Returns length of thing */
  int i=0, c;

  while ( (c=fgetc(file)) != SEPARATOR )
    if (c==EOF)
      return 0;
    else
      string[i++]=c;
  string[i]='\0';
  return i;
} /* fgetz */

static void skip_comment(FILE *file) {
  int i;
  
  do  { 
    i=getc(file);

    if(i=='#')
      do 
	i=getc(file);
      while(i != '\n');
  } while(isspace(i));

  ungetc(i, file);

  return;
} /* skip_comment */
