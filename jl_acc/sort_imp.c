
/* this file possibly included multiple times by sort.c
   for sorting different integer sizes;
   
   look in sort.c for some controlling macro definitions,
   like Value, Index, Data, and function names */

#ifdef Value

/*------------------------------------------------------------------------------
  
  Radix Sort
  
  stable; O(n+k) time and extra storage
    where k = (digits in an int) * 2^(bits per digit)
    (e.g. k = 4 * 256 = 1024 for 32-bit ints with 8-bit digits)

  brief description:
    input sorted stably on each digit, starting with the least significant
    counting sort is used for each digit:
      a pass through the input counts the occurences of each digit value
      on a second pass, each input has a known destination
  
  tricks:
    all counting passes are combined into one
    the counting pass also computes the inclusive bit-wise or of all inputs,
      which is used to skip digit positions for which all inputs have zeros

  ----------------------------------------------------------------------------*/

#define DIGIT_BITS   8
#define DIGIT_VALUES (1<<DIGIT_BITS)
#define DIGIT_MASK   ((Value)(DIGIT_VALUES-1))
#define CEILDIV(a,b) (((a)+(b)-1)/(b))
#define DIGITS       CEILDIV(CHAR_BIT*sizeof(Value),DIGIT_BITS)
#define VALUE_BITS   (DIGIT_BITS*DIGITS)
#define COUNT_SIZE   (DIGITS*DIGIT_VALUES)

/* used to unroll a tiny loop: */
#define COUNT_DIGIT_01(n,i) \
    if(n>i) count[i][val&DIGIT_MASK]++, val>>=DIGIT_BITS
#define COUNT_DIGIT_02(n,i) COUNT_DIGIT_01(n,i); COUNT_DIGIT_01(n,i+ 1)
#define COUNT_DIGIT_04(n,i) COUNT_DIGIT_02(n,i); COUNT_DIGIT_02(n,i+ 2)
#define COUNT_DIGIT_08(n,i) COUNT_DIGIT_04(n,i); COUNT_DIGIT_04(n,i+ 4)
#define COUNT_DIGIT_16(n,i) COUNT_DIGIT_08(n,i); COUNT_DIGIT_08(n,i+ 8)
#define COUNT_DIGIT_32(n,i) COUNT_DIGIT_16(n,i); COUNT_DIGIT_16(n,i+16)
#define COUNT_DIGIT_64(n,i) COUNT_DIGIT_32(n,i); COUNT_DIGIT_32(n,i+32)

static Value radix_count(const Value *A, const Value *end, Index stride,
                         Index count[DIGITS][DIGIT_VALUES])
{
  Value bitorkey = 0;
  memset(count,0,COUNT_SIZE*sizeof(Index));
  do {
    Value val=*A;
    bitorkey|=val;
    COUNT_DIGIT_64(DIGITS,0);
    /* above macro expands to:
    if(DIGITS> 0) count[ 0][val&DIGIT_MASK]++, val>>=DIGIT_BITS;
    if(DIGITS> 1) count[ 1][val&DIGIT_MASK]++, val>>=DIGIT_BITS;
      ...
    if(DIGITS>63) count[63][val&DIGIT_MASK]++, val>>=DIGIT_BITS;
    */
  } while(A+=stride,A!=end);
  return bitorkey;
}

#undef COUNT_DIGIT_01
#undef COUNT_DIGIT_02
#undef COUNT_DIGIT_04
#undef COUNT_DIGIT_08
#undef COUNT_DIGIT_16
#undef COUNT_DIGIT_32
#undef COUNT_DIGIT_64

static void radix_offsets(Index *c)
{
  Index sum=0, t, *ce=c+DIGIT_VALUES;
  do t=*c, *c++ = sum, sum+=t; while(c!=ce);
}

static unsigned radix_zeros(Value bitorkey, Index count[DIGITS][DIGIT_VALUES],
                            unsigned *shift, Index **offsets)
{
  unsigned digits=0, sh=0; Index *c = &count[0][0];
  do {
    if(bitorkey&DIGIT_MASK) *shift++ = sh, *offsets++ = c, ++digits,
                            radix_offsets(c);
  } while(bitorkey>>=DIGIT_BITS,sh+=DIGIT_BITS,c+=DIGIT_VALUES,sh!=VALUE_BITS);
  return digits;
}

static void radix_pass(const Value *A, const Value *end, Index stride,
                       unsigned sh, Index *off, Value *out)
{
  do out[off[(*A>>sh)&DIGIT_MASK]++] = *A; while(A+=stride,A!=end);
}

static void radix_sort(const Value *A, Index n, Index stride,
                       Value *out, Value *work)
{
  Index count[DIGITS][DIGIT_VALUES];
  const Value *end = A+n*stride;
  Value bitorkey = radix_count(A, end, stride, count);
  unsigned shift[DIGITS]; Index *offsets[DIGITS];
  unsigned digits = radix_zeros(bitorkey,count,shift,offsets);
  if(digits==0) {
    memset(out,0,sizeof(Value)*n);
  } else {
    Value *src, *dst; unsigned d;
    if((digits&1)==1) src=out,dst=work;
                 else dst=out,src=work;
    radix_pass(A,end,stride,shift[0],offsets[0],src);
    for(d=1;d!=digits;++d) {
      Value *t;
      radix_pass(src,src+n,1,shift[d],offsets[d],dst);
      t=src,src=dst,dst=t;
    }
  }
}

static void radix_index_pass_b(const Value *A, Index n, Index stride,
                               unsigned sh, Index *off, Data *out)
{
  Index i=0;
  do {
    Value v = *A;
    Data *d = &out[off[(v>>sh)&DIGIT_MASK]++];
    d->v=v, d->i=i++;
  } while(A+=stride,i!=n);
}

static void radix_index_pass_m(const Data *src, const Data *end,
                               unsigned sh, Index *off, Data *out)
{
  do {
    Data *d = &out[off[(src->v>>sh)&DIGIT_MASK]++];
    d->v=src->v,d->i=src->i;
  } while(++src!=end);
}

static void radix_index_pass_e(const Data *src, const Data *end,
                               unsigned sh, Index *off,
                               Index *out)
{
  do out[off[(src->v>>sh)&DIGIT_MASK]++]=src->i; while(++src!=end);
}

static void radix_index_pass_be(const Value *A, Index n, Index stride,
                                unsigned sh, Index *off, Index *out)
{
  Index i=0;
  do out[off[(*A>>sh)&DIGIT_MASK]++]=i++; while(A+=stride,i!=n);
}

static void radix_index_sort(const Value *A, Index n, Index stride,
                             Index *idx, Data *work)
{
  Index count[DIGITS][DIGIT_VALUES];
  Value bitorkey = radix_count(A, A+n*stride, stride, count);
  unsigned shift[DIGITS]; Index *offsets[DIGITS];
  unsigned digits = radix_zeros(bitorkey,count,shift,offsets);
  if(digits==0) {
    Index i=0; do *idx++=i++; while(i!=n);
  } else if(digits==1) {
    radix_index_pass_be(A,n,stride,shift[0],offsets[0],idx);
  } else {
    Data *src, *dst; unsigned d;
    if((digits&1)==0) dst=work,src=dst+n;
                 else src=work,dst=src+n;
    radix_index_pass_b(A,n,stride,shift[0],offsets[0],src);
    for(d=1;d!=digits-1;++d) {
      Data *t;
      radix_index_pass_m(src,src+n,shift[d],offsets[d],dst);
      t=src,src=dst,dst=t;
    }
    radix_index_pass_e(src,src+n,shift[d],offsets[d],idx);
  }
}

/*------------------------------------------------------------------------------
  
  Merge Sort
  
  stable; O(n log n) time

  ----------------------------------------------------------------------------*/

static void merge_sort(const Value *A, Index n, Index stride,
                       Value *out, Value *work)
{
  Value *buf[2]={out,work};
  Index base=-n, odd=0, c=0, b=1;
  for(;;) {
    Value *p;
    if((c&1)==0) {
      base+=n, n+=(odd&1), c|=1, b^=1;
      while(n>3) odd<<=1,odd|=(n&1),n>>=1,c<<=1,b^=1;
    } else
      base-=n-(odd&1),n<<=1,n-=(odd&1),odd>>=1,c>>=1;
    if(c==0) break;
    p = buf[b]+base;
    if(n==2) {
      Value v[2]; v[0]=*A,A+=stride,v[1]=*A,A+=stride;
      if(v[1]<v[0]) p[0]=v[1],p[1]=v[0];
               else p[0]=v[0],p[1]=v[1];
    } else if(n==3) {
      Value v[3]; v[0]=*A,A+=stride,v[1]=*A,A+=stride,v[2]=*A,A+=stride;
      if(v[1]<v[0]) {
        if(v[2]<v[1])        p[0]=v[2],p[1]=v[1],p[2]=v[0];
        else { if(v[2]<v[0]) p[0]=v[1],p[1]=v[2],p[2]=v[0];
                        else p[0]=v[1],p[1]=v[0],p[2]=v[2]; }
      } else {
        if(v[2]<v[0])        p[0]=v[2],p[1]=v[0],p[2]=v[1];
        else { if(v[2]<v[1]) p[0]=v[0],p[1]=v[2],p[2]=v[1];
                        else p[0]=v[0],p[1]=v[1],p[2]=v[2]; }
      }
    } else {
      const Index na = n>>1, nb = (n+1)>>1;
      const Value *ap = buf[b^1]+base, *ae = ap+na;
      Value *bp = p+na, *be = bp+nb;
      for(;;) {
        if(*bp<*ap) {
          *p++=*bp++;
          if(bp!=be) continue;
          do *p++=*ap++; while(ap!=ae);
          break;
        } else {
          *p++=*ap++;
          if(ap==ae) break;
        }
      }
    }
  }
}

static void merge_index_sort(const Value *A, const Index An, Index stride,
                             Index *idx, Data *work)
{
  Data *buf[2]={work+An,work};
  Index n=An, base=-n, odd=0, c=0, b=1;
  Index i=0;
  for(;;) {
    Data *p;
    if((c&1)==0) {
      base+=n, n+=(odd&1), c|=1, b^=1;
      while(n>3) odd<<=1,odd|=(n&1),n>>=1,c<<=1,b^=1;
    } else
      base-=n-(odd&1),n<<=1,n-=(odd&1),odd>>=1,c>>=1;
    if(c==0) break;
    p = buf[b]+base;
    if(n==2) {
      Value v[2]; v[0]=*A,A+=stride,v[1]=*A,A+=stride;
      if(v[1]<v[0]) p[0].v=v[1],p[0].i=i+1, p[1].v=v[0],p[1].i=i  ;
               else p[0].v=v[0],p[0].i=i  , p[1].v=v[1],p[1].i=i+1;
      i+=2;
    } else if(n==3) {
      Value v[3]; v[0]=*A,A+=stride,v[1]=*A,A+=stride,v[2]=*A,A+=stride;
      if(v[1]<v[0]) {
        if(v[2]<v[1])        p[0].v=v[2],p[1].v=v[1],p[2].v=v[0],
                             p[0].i=i+2 ,p[1].i=i+1 ,p[2].i=i   ;
        else { if(v[2]<v[0]) p[0].v=v[1],p[1].v=v[2],p[2].v=v[0],
                             p[0].i=i+1 ,p[1].i=i+2 ,p[2].i=i   ;
                        else p[0].v=v[1],p[1].v=v[0],p[2].v=v[2],
                             p[0].i=i+1 ,p[1].i=i   ,p[2].i=i+2 ; }
      } else {
        if(v[2]<v[0])        p[0].v=v[2],p[1].v=v[0],p[2].v=v[1],
                             p[0].i=i+2 ,p[1].i=i   ,p[2].i=i+1 ;
        else { if(v[2]<v[1]) p[0].v=v[0],p[1].v=v[2],p[2].v=v[1],
                             p[0].i=i   ,p[1].i=i+2 ,p[2].i=i+1 ;
                        else p[0].v=v[0],p[1].v=v[1],p[2].v=v[2],
                             p[0].i=i   ,p[1].i=i+1 ,p[2].i=i+2 ; }
      }
      i+=3;
    } else {
      const Index na = n>>1, nb = (n+1)>>1;
      const Data *ap = buf[b^1]+base, *ae = ap+na;
      Data *bp = p+na, *be = bp+nb;
      for(;;) {
        if(bp->v<ap->v) {
          *p++=*bp++;
          if(bp!=be) continue;
          do *p++=*ap++; while(ap!=ae);
          break;
        } else {
          *p++=*ap++;
          if(ap==ae) break;
        }
      }
    }
  }
  {
    const Data *p = buf[0], *pe = p+An;
    do *idx++ = (p++)->i; while(p!=pe);
  }
}

/*------------------------------------------------------------------------------
  
  Hybrid Stable Sort
  
  low-overhead merge sort when n is small,
  otherwise asymptotically superior radix sort

  result = O(n) sort with good performance for all n
  
  A, n, stride : specifices the input
  
  sort:
     Value out[n] : the sorted values (output)
     Value work[n]: scratch area
  
  index_sort:
     Index idx[n]  : the sorted indices (output)
     Data work[2*n]: scratch area

  ----------------------------------------------------------------------------*/

void sort(const Value *A, Index n, Index stride, Value *out, Value *work)
{
  if(n<DIGIT_VALUES) {
    if(n==0) return;
    if(n==1) *out = *A;
    else     merge_sort(A,n,stride,out,work);
  } else     radix_sort(A,n,stride,out,work);
}

void index_sort(const Value *A, Index n, Index stride,
                Index *idx, Data *work)
{
  if(n<DIGIT_VALUES) {
    if(n==0) return;
    if(n==1) *idx=0;
    else     merge_index_sort(A,n,stride,idx,work);
  } else     radix_index_sort(A,n,stride,idx,work);
}

#undef DIGIT_BITS
#undef DIGIT_VALUES
#undef DIGIT_MASK
#undef CEILDIV
#undef DIGITS
#undef VALUE_BITS
#undef COUNT_SIZE

#endif
