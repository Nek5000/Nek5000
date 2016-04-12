#include <ctype.h>
#include "iniparser.h"

#ifndef FNAME_H
#define FNAME_H

/*
   FORTRAN naming convention
     default      cpgs_setup, etc.
     -DUPCASE     CPGS_SETUP, etc.
     -DUNDERSCORE cpgs_setup_, etc.
*/

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#endif

#define finiparser_dump         FORTRAN_NAME(finiparser_dump,  FINIPARSER_DUMP)
#define finiparser_findInvalid  FORTRAN_NAME(finiparser_findinvalid,  FINIPARSER_FINDINVALID)
#define finiparser_load         FORTRAN_NAME(finiparser_load,  FINIPARSER_LOAD)
#define finiparser_free         FORTRAN_NAME(finiparser_free,  FINIPARSER_FREE)
#define finiparser_getString    FORTRAN_NAME(finiparser_getstring,  FINIPARSER_GETSTRING)
#define finiparser_getBool      FORTRAN_NAME(finiparser_getbool,  FINIPARSER_GETBOOL)
#define finiparser_find         FORTRAN_NAME(finiparser_find,  FINIPARSER_FIND)
#define finiparser_getDbl       FORTRAN_NAME(finiparser_getdbl,  FINIPARSER_GETDBL)
#define finiparser_getToken     FORTRAN_NAME(finiparser_gettoken,  FINIPARSER_GETTOKEN)
#define finiparser_findTokens   FORTRAN_NAME(finiparser_findtokens,  FINIPARSER_FINDTOKENS)

#define ntokenmax  100 

static dictionary *dic=NULL;
static char *token[ntokenmax];


char *addchar0(char * str,int str_len)
{
    int i, real_len;
    char *newstr;

    /* strip trailing blanks in datarep */
    if (str <= (char *) 0) {
        return NULL;
    }
    for (i=str_len-1; i>=0; i--) if (str[i] != ' ') break;
    if (i < 0) {
        return NULL;
    }
    real_len = i + 1;

    newstr = (char *) malloc((real_len+1)*sizeof(char));
    strncpy(newstr, str, real_len);
    newstr[real_len] = '\0';
    return newstr;
}

void finiparser_dump()
{
    if(dic != NULL) iniparser_dump(dic,stdout);
    return;
}

void finiparser_find(int* out,char *key,int* ifnd,int key_len)
{
    int tmp;
    *ifnd = 0;
    tmp = iniparser_find_entry(dic,addchar0(key,key_len));
    if (tmp == 1) {
       *out = tmp;
       *ifnd = 1;
    }
    return;
}

void finiparser_findInvalid(int* ifnd)
{
    int i;

    *ifnd = 0;

/* cannot be done w/o a reference dictionary containing all the valid keys
    for (i=0 ; i < dic->size ; i++) {
        if (dic->key[i]==NULL) continue ;
        printf("%s = [%s]\n", dic->key[i], dic->val[i]);
        if (iniparser_find_entry(dic_ref,dic->key[i]) == 0) {
            printf(" Invalid entry %s = [%s]\n", dic_->key[i], dic->val[i]);
            *ifnd++;
        } 
    }
*/
    return;
}

void finiparser_load(char * fname,int* ierr,int fname_len)
{
    *ierr = 0;
    dic = iniparser_load(addchar0(fname,fname_len)); 
    if (dic == NULL) *ierr = 1;
    return;
}

void finiparser_free()
{
    dictionary_del(dic);
    return;
}

void finiparser_getString(char *out,char *key,int *ifnd,int out_len,int key_len)
{
    int i;
    const char* str;
    int real_out_len;

    *ifnd = 0;
    for (i=0; i<out_len; i++) out[i] = ' ';

    str = iniparser_getstring(dic,addchar0(key,key_len),NULL);
    if (str != NULL) {
       real_out_len = strlen(str);
       if(real_out_len <= out_len) {
         strncpy(out,str,real_out_len);
         *ifnd = 1;
       } 
    }
    return;
}

void finiparser_getBool(int* out,char *key,int* ifnd,int key_len)
{
    int tmp;

    *ifnd = 0;
    tmp = iniparser_getboolean(dic,addchar0(key,key_len),-1);
    if (tmp != -1) { 
       *out = tmp;
       *ifnd = 1;
    }
    return;
}

void finiparser_getDbl(double* out,char *key,int *ifnd,int key_len)
{
    const char* str;

    *ifnd = 0;
    str = iniparser_getstring(dic,addchar0(key,key_len),NULL);
    if (str != NULL) {
       *out = atof(str);
       *ifnd = 1;
    }
    return;
}

void finiparser_getToken(char *out,int *id,int out_len)
{
    int real_out_len, i;

    for (i=0; i<out_len; i++) out[i] = ' ';
    if(*id > ntokenmax) return;
    real_out_len = strlen(token[*id-1]);
    if(real_out_len <= out_len) strncpy(out,token[*id-1],real_out_len);
    return;
}

void finiparser_findTokens(char *key, char *delim, int *icounter,int key_len,int delim_len)
{
    char *str;
    char *d;
    int i = 0;

    *icounter = 0;

    d = addchar0(delim,delim_len);
    str = iniparser_getstring(dic,addchar0(key,key_len),NULL);
    if (str == NULL) return;

    token[0] = strtok(str,d);
    while (token[i] != NULL && i <= ntokenmax-1) token[++i] = strtok(NULL,d);
    *icounter = i;
    
   return;
}
