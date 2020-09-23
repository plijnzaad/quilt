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

/* module for doing simple hash table stuff. Both keys and values are
 * either fixed length, or null terminated strings; this is specified
 * when setting up the table, with Hcreate. A key_size of zero implies a
 * 0-terminated string. A key_size of -1 signifies a pointer that is to be
 * hashed itself, instead of its contents (for quick comparison of pointers)
 */

#include <stdlib.h>			/* for definition of int */

#ifndef _HASH_H_
#define _HASH_H_
#include "compat.h"

typedef struct hnode_ {			/* holds a key/value pair */
  struct hnode_ * next;			/* make it resemble a 'list_p' */
  void * key;
  void * value;
  struct hnode_ * nexthnode;		/* for internal use */

#ifdef FASTHASH
  unsigned int second_hash;		/* saves calls to strcmp/memcmp */
#endif

} hnode_t, *hnode_p;

typedef struct htable_ {		/* user: only read contents */
  int size;				/* size of table (nr of slots) */
  hnode_p * table;			/* the table */
  int key_size;				/* length in bytes of keys; if 0 */
/*  int valuesize; */			/* '\0' terminated strings */
  int (*difference_test)();		/* either strcmp or memcmp */
  int nkeys;				/* total nr of keys */
  int unused;				/* nr of unused slots */

  /* following stuff is only valid if keys != NULL (after a call to Hkeys) */
  hnode_p keys;				/* NULL if stale; otherwise list */
  hnode_p sorted_keys;			/* NULL if stale; otherwise list */
  int maxlen;				/* max list length of one slot */
} htable_t, *htable_p;

htable_p Hcreate(htable_p init, int tablesize, int key_size);
/* sets up htable of size (actually, the first prime above it). Loads it with
 * init (unless NULL, and it must be a hashtable with the same keysize).
 * Loading does a Hfind on all the table entries from init, and copies the
 * value pointers (so both tables refer to the same value data. To obtain
 * indivual data, consider 
 *   for(node=Hkeys(ht); node; node=b->next)  
 *     node->value=memsav(node->value, value_size);
 * A pointer to the newly created hashtable is returned
 */

#define Hfree(T, F) _Hfree((T), (F), __FILE__, __LINE__)
void _Hfree(htable_p, void(*free_func)(), const char *, int);
  /* deletes complete htable, using free_func (unless NULL) on hnode->value
   * If the key_size < 0, the free_func (if any) is called with the KEY as
   * argument, freeing a pointer that was not allocated in hash.c ! (this is
   * intended to be so) 
   */

hnode_p Hget(htable_p table, const void *key);
/* sees if KEY is defined in TABLE; returns NULL if not found. Otherwise
 * returns pointer to entry. Does NOT create an entry if it's not
 * found. To speed up searching, a found entry is moved to the head of the
 * list forming the bucket, so the ordering of a Hkeys() call is lost
 */

hnode_p Hfeel(htable_p, const void *key);
/* like Hgetined, but doesn't reorder, so leaves the ordering of a Hkeys intact
 */

hnode_p _Hset(htable_p, const void *key, const char * file, int linenr);
#define Hset(T, K) _Hset((T), (K), __FILE__, __LINE__)
/* searches key in htable; creates one if not found. Returns pointer to
 * entry. After creation of the key, the contents of the returned pointer
 * are garantueed to be NULL.
 */

void * _Hdelete(htable_p, const void *key, const char *, int);
#define Hdelete(T,K) _Hdelete((T),(K), __FILE__, __LINE__)
  /* deletes key from htable. Returns the value pointer, which may need */
  /* freeing itself */

hnode_p Hkeys(htable_p);
  /* returns keys of htable_p in apparently random order (but it can also */
  /* be the ordering from a previous sort). The returned list becomes */
  /* invalid after calls to Hgetined, Hfind, Hdelete; after Hsorted_keys, */
  /* the ordering may be different. After Hfeel, the order is the same */

hnode_p Hsorted_keys(htable_p, int (*compare)(void*, void*));
  /* returns keys of htable_p in order sorted by compare function. If */
  /* compare == NULL, comparison will be done by strcmp on the keys for */
  /* key_size==0, by memcmp otherwise. and by 2*((a < b)-1)*(a!=b) if -1 */
  /* The compare function, with type 'int (*compare)(hnode_p * one, */
  /* hnode_p * other)', compares pointers to hnode_p's, not hnode_p's */
  /* themselves !. The returned list becomes invalid after calls to */
  /* Hget, Hset, and Hdelete */ 

int Hread(FILE*file, htable_p table, char *keyfmt, char *valfmt, 
	   int valsize, 
	  void*(*parse_func)(const char *string));
/* reads data from FILE into TABLE, under control of KEYFMT and VALFMT,
 * and assuming a value size of VALSIZE bytes. If KEYFMT is not null,
 * comments (# ... \n) and empty lines are skipped. Otherwise, items are
 * read as appropriate to their size, or, if size is zero, up to the next
 * separating character (SEPARTOR in hash.c, usually just '\0'). If
 * PARSE_FUNC is NULL, all node->value's are saved with memsav(value).
 * Otherwise they are set to PARSE_FUNC(value); ie. PARSE_FUNC parses and
 * allocates the values to be kept. (with return value NULL implying an
 * error) Hread returns number of new entries.
 */

void Hwrite(FILE*file, htable_p table, char *keyfmt, char *valfmt,
	    int valsize);
/* writes table to file, under control of keyfmt and valfmt, and
 * assuming a value size of valsize bytes. If an item's fmt is NULL or "", AND 
 * its size is 0, a '\0' is used as separator character.
 */

#endif /* !defined _HASH_H_ */
