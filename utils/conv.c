#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h> 

#define ERROR do {fprintf(stderr, \
			"Usage: conv [0x<n> | 0<n> | <n> | \'<c>\']\n");\
		  return 1; } while (0)

#define BIT(N) (1<<(N))			/* BIT(0)=1, BIT(1)=2, BIT(3)=4,etc. */

#define ZERO_CHAR '0'			/* used to be '.'  */
#define ONE_CHAR '1'			/* used to be '|' */

#include <limits.h> 
#if defined LONG_BIT
#	define SIZE LONG_BIT
#elif defined WORD_BIT			/* determine word length */
#	define SIZE WORD_BIT
#else
#	define SIZE 32
#endif

char * intbits( unsigned long I) {
  int i; 
  char *cp;
  static char word[ SIZE+1 ];

  word[ SIZE ]='\0';			/* terminating 0 */
  for (i=0,cp = &(word[ SIZE-1 ]); i < SIZE ; i++, cp--)
    *cp = (I & BIT(i) )? ONE_CHAR : ZERO_CHAR ;
  return cp+1;				/* since decremented one too far */
} /* intbits */

static int find_set_bits(char *setbits, const char *bits)  {
  /* writes the offsets of bits that are set into intbits; returns total nr */
  int i, n;
  const char *cp;

  i=n=0;
  cp= bits + strlen(bits)-1;		/* start at end */
  while(cp >= bits) { 
    if( *cp == ONE_CHAR )
      sprintf(setbits, "%s %d ", setbits, i), n++;
    cp--, i++;
  }
  return n;
} /* find_set_bits */

int convert(char * string) {
  long i;
  int n;
  char * cp;
  char bits[80], setbits[240];

  bzero(setbits, sizeof(setbits));

  i=strtoul(string, &cp, 0);
  if (*cp) {				/* maybe something wrong */
    if(cp==string) {			/* i.e., not recognized at all */
      /* it still may be one letter, the ASCII value of which is desired */
      if ( isspace(cp[1]) || cp[1]==0 )	/* may yet be single character */
	i=(int)cp[0];
      else
	ERROR;
    } else 
      ERROR;
  }
  strcpy(bits, intbits(i));
  n=find_set_bits(setbits, bits);
  printf("%ld (%lu) 0%lo 0x%lx ; %d bits set: %s%s\n\t%s ``%c\'\'\n",
	   i, (unsigned long)i, i, i, 
	 n, (n<=6)?"":"\n", setbits, bits, (char)i);

  return 0;
} /* convert */

/* outputs number as decimal, octal and hexadecimal, bitpattern and char */
int main(int argc, char ** argv) {
  int err=0;
  char line[ BUFSIZ ];			/* lines from stdin */

  if (argc != 1) {			/* a couple of args */
    argv++;
    while(--argc)
      err += convert(*argv++);
    return err;
  }

  while(gets(line))			/* else: take lines from stdin */
    err += convert(line);
  return err;
} /* main */
