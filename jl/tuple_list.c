#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "errmem.h"
#include "types.h"
#include "minmax.h"
#include "sort.h"

typedef struct {
  unsigned mi,ml,mr;
  uint n, max;
  sint *vi; slong *vl; real *vr;
} tuple_list;

void tuple_list_permute(tuple_list *tl, uint *perm, void *work)
{
  const unsigned mi=tl->mi, ml=tl->ml, mr=tl->mr;
  const unsigned int_size  = mi*sizeof(sint),
                 long_size = ml*sizeof(slong),
                 real_size = mr*sizeof(real);
  if(mi) {
    uint *p=perm, *pe=p+tl->n; char *sorted=work;
    while(p!=pe) memcpy(sorted,&tl->vi[mi*(*p++)],int_size),sorted+=int_size;
    memcpy(tl->vi,work,int_size*tl->n);
  }
  if(ml) {
    uint *p=perm, *pe=p+tl->n; char *sorted=work;
    while(p!=pe) memcpy(sorted,&tl->vl[ml*(*p++)],long_size),sorted+=long_size;
    memcpy(tl->vl,work,long_size*tl->n);
  }
  if(mr) {
    uint *p=perm, *pe=p+tl->n; char *sorted=work;
    while(p!=pe) memcpy(sorted,&tl->vr[mr*(*p++)],real_size),sorted+=real_size;
    memcpy(tl->vr,work,real_size*tl->n);
  }
}

void tuple_list_sort(tuple_list *tl, unsigned key, buffer *buf)
{
  const unsigned mi=tl->mi, ml=tl->ml, mr=tl->mr;
  const unsigned int_size =  mi*sizeof(sint);
  const unsigned long_size = ml*sizeof(slong);
  const unsigned real_size = mr*sizeof(real);
  const unsigned width = umax_3(int_size,long_size,real_size);
  const unsigned data_size = key>=mi ? sizeof(sort_data_long):sizeof(sort_data);
  uint work_min=tl->n * umax_2(2*data_size,sizeof(sint)+width);
  uint *work;
  buffer_reserve(buf,work_min);
  work = buf->ptr;
  if(key<mi)
    index_sort     ((uint *)&tl->vi[key   ],tl->n,mi, work, (void*)work);
  else
    index_sort_long((ulong*)&tl->vl[key-mi],tl->n,ml, work, (void*)work);
  tuple_list_permute(tl,work,work+tl->n);
}

