#if !defined(T) || !defined(SORT_SUFFIX)
#error sort_imp.c not meant to be compiled by itself
#endif

#define sort_data       TOKEN_PASTE(sort_data      ,SORT_SUFFIX)
#define radix_count     TOKEN_PASTE(radix_count    ,SORT_SUFFIX)
#define radix_offsets   TOKEN_PASTE(radix_offsets  ,SORT_SUFFIX)
#define radix_zeros     TOKEN_PASTE(radix_zeros    ,SORT_SUFFIX)
#define radix_passv     TOKEN_PASTE(radix_passv    ,SORT_SUFFIX)
#define radix_sortv     TOKEN_PASTE(radix_sortv    ,SORT_SUFFIX)
#define radix_passp0_b  TOKEN_PASTE(radix_passp0_b ,SORT_SUFFIX)
#define radix_passp_b   TOKEN_PASTE(radix_passp_b  ,SORT_SUFFIX)
#define radix_passp_m   TOKEN_PASTE(radix_passp_m  ,SORT_SUFFIX)
#define radix_passp_e   TOKEN_PASTE(radix_passp_e  ,SORT_SUFFIX)
#define radix_passp0_be TOKEN_PASTE(radix_passp0_be,SORT_SUFFIX)
#define radix_passp_be  TOKEN_PASTE(radix_passp_be, SORT_SUFFIX)
#define radix_sortp     TOKEN_PASTE(radix_sortp    ,SORT_SUFFIX)
#define merge_sortv     TOKEN_PASTE(merge_sortv    ,SORT_SUFFIX)
#define merge_sortp0    TOKEN_PASTE(merge_sortp0   ,SORT_SUFFIX)
#define merge_sortp     TOKEN_PASTE(merge_sortp    ,SORT_SUFFIX)

#ifdef PREFIX
#  define sortv TOKEN_PASTE(TOKEN_PASTE(PREFIX,sortv),SORT_SUFFIX)
#  define sortp TOKEN_PASTE(TOKEN_PASTE(PREFIX,sortp),SORT_SUFFIX)
#else
#  define sortv TOKEN_PASTE(sortv,SORT_SUFFIX)
#  define sortp TOKEN_PASTE(sortp,SORT_SUFFIX)
#endif

typedef struct { T v; uint i; } sort_data;

#define INC_PTR(A,stride) ((A)=(T*)((char*)(A)+(stride)))
#define INDEX_PTR(A,stride,i) (*(T*)((char*)(A)+(i)*(stride)))

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
#define DIGIT_MASK   ((T)(DIGIT_VALUES-1))
#define CEILDIV(a,b) (((a)+(b)-1)/(b))
#define DIGITS       CEILDIV(CHAR_BIT*sizeof(T),DIGIT_BITS)
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

static T radix_count(uint count[DIGITS][DIGIT_VALUES],
                     const T *A, const T *end, const unsigned stride)
{
  T bitorkey = 0;
  memset(count,0,COUNT_SIZE*sizeof(uint));
  do {
    T val=*A;
    bitorkey|=val;
    COUNT_DIGIT_64(DIGITS,0);
    /* above macro expands to:
    if(DIGITS> 0) count[ 0][val&DIGIT_MASK]++, val>>=DIGIT_BITS;
    if(DIGITS> 1) count[ 1][val&DIGIT_MASK]++, val>>=DIGIT_BITS;
      ...
    if(DIGITS>63) count[63][val&DIGIT_MASK]++, val>>=DIGIT_BITS;
    */
  } while(INC_PTR(A,stride),A!=end);
  return bitorkey;
}

#undef COUNT_DIGIT_01
#undef COUNT_DIGIT_02
#undef COUNT_DIGIT_04
#undef COUNT_DIGIT_08
#undef COUNT_DIGIT_16
#undef COUNT_DIGIT_32
#undef COUNT_DIGIT_64

static void radix_offsets(uint *c)
{
  uint sum=0, t, *ce=c+DIGIT_VALUES;
  do t=*c, *c++ = sum, sum+=t; while(c!=ce);
}

static unsigned radix_zeros(T bitorkey, uint count[DIGITS][DIGIT_VALUES],
                            unsigned *shift, uint **offsets)
{
  unsigned digits=0, sh=0; uint *c = &count[0][0];
  do {
    if(bitorkey&DIGIT_MASK) *shift++ = sh, *offsets++ = c, ++digits,
                            radix_offsets(c);
  } while(bitorkey>>=DIGIT_BITS,sh+=DIGIT_BITS,c+=DIGIT_VALUES,sh!=VALUE_BITS);
  return digits;
}

static void radix_passv(const T *A, const T *end, const unsigned stride,
                        const unsigned sh, uint *off, T *const out)
{
  do out[off[(*A>>sh)&DIGIT_MASK]++] = *A; while(INC_PTR(A,stride),A!=end);
}

static void radix_sortv(T *out, const T *A, uint n, const unsigned stride,
                        T *work)
{
  uint count[DIGITS][DIGIT_VALUES];
  const T *end = &INDEX_PTR(A,stride,n);
  T bitorkey = radix_count(count, A,end,stride);
  unsigned shift[DIGITS]; uint *offsets[DIGITS];
  unsigned digits = radix_zeros(bitorkey,count,shift,offsets);
  if(digits==0) {
    memset(out,0,n*sizeof(T));
  } else {
    T *src, *dst; unsigned d;
    if((digits&1)==1) src=out,dst=work;
                 else dst=out,src=work;
    radix_passv(A,end,stride,shift[0],offsets[0],src);
    for(d=1;d!=digits;++d) {
      T *t;
      radix_passv(src,src+n,sizeof(T),shift[d],offsets[d],dst);
      t=src,src=dst,dst=t;
    }
  }
}

static void radix_passp0_b(const T *A, uint n, const unsigned stride,
                           const unsigned sh, uint *const off,
                           sort_data *const out)
{
  uint i=0;
  do {
    T v = *A;
    sort_data *d = &out[off[(v>>sh)&DIGIT_MASK]++];
    d->v=v, d->i=i++;
  } while(INC_PTR(A,stride),i!=n);
}

static void radix_passp_b(const uint *p,
                          const T *A, uint n, const unsigned stride,
                          const unsigned sh, uint *const off,
                          sort_data *const out)
{
  const uint *pe = p+n;
  do {
    uint j = *p++;
    T v = INDEX_PTR(A,stride,j);
    sort_data *d = &out[off[(v>>sh)&DIGIT_MASK]++];
    d->v=v, d->i=j;
  } while(p!=pe);
}

static void radix_passp_m(const sort_data *src, const sort_data *end,
                          const unsigned sh, uint *const off,
                          sort_data *const out)
{
  do {
    sort_data *d = &out[off[(src->v>>sh)&DIGIT_MASK]++];
    d->v=src->v,d->i=src->i;
  } while(++src!=end);
}

static void radix_passp_e(const sort_data *src, const sort_data *end,
                          const unsigned sh, uint *const off,
                          uint *const out)
{
  do out[off[(src->v>>sh)&DIGIT_MASK]++]=src->i; while(++src!=end);
}

static void radix_passp0_be(uint *const out,
                            const T *A, uint n, const unsigned stride,
                            const unsigned sh, uint *const off)
{
  uint i=0;
  do out[off[(*A>>sh)&DIGIT_MASK]++]=i++; while(INC_PTR(A,stride),i!=n);
}

static void radix_passp_be(uint *p,
                           const T *A, uint n, const unsigned stride,
                           const unsigned sh, uint *const off,
                           sort_data *work)
{
  uint *q = p, *qe = p+n;
  do {
    uint j = *q++;
    T v = INDEX_PTR(A,stride,j);
    work[off[(v>>sh)&DIGIT_MASK]++].i=j;
  } while(q!=qe);
  q=p,qe=p+n; do *q++ = (work++)->i; while(q!=qe);
}

static void radix_sortp(uint *idx, uint perm_start,
                        const T *A, uint n, const unsigned stride,
                        sort_data *work)
{
  uint count[DIGITS][DIGIT_VALUES];
  T bitorkey = radix_count(count, A,&INDEX_PTR(A,stride,n),stride);
  unsigned shift[DIGITS]; uint *offsets[DIGITS];
  unsigned digits = radix_zeros(bitorkey,count,shift,offsets);
  if(digits==0) {
    if(!perm_start) { uint i=0; do *idx++=i++; while(i!=n); }
  } else if(digits==1) {
    if(perm_start) radix_passp_be (idx,A,n,stride,shift[0],offsets[0],work);
              else radix_passp0_be(idx,A,n,stride,shift[0],offsets[0]);
  } else {
    sort_data *src, *dst; unsigned d;
    if((digits&1)==0) dst=work,src=dst+n;
                 else src=work,dst=src+n;
    if(perm_start) radix_passp_b (idx,A,n,stride,shift[0],offsets[0],src);
              else radix_passp0_b(    A,n,stride,shift[0],offsets[0],src);
    for(d=1;d!=digits-1;++d) {
      sort_data *t;
      radix_passp_m(src,src+n,shift[d],offsets[d],dst);
      t=src,src=dst,dst=t;
    }
    radix_passp_e(src,src+n,shift[d],offsets[d],idx);
  }
}

/*------------------------------------------------------------------------------
  
  Merge Sort
  
  stable; O(n log n) time

  ----------------------------------------------------------------------------*/

#define MERGE_2(p,v)                           \
  if(VAL(v[1])<VAL(v[0])) p[0]=v[1],p[1]=v[0]; \
                     else p[0]=v[0],p[1]=v[1]
#define MERGE_3(p,v) do                                              \
  if(VAL(v[1])<VAL(v[0])) {                                          \
    if(VAL(v[2])<VAL(v[1]))        p[0]=v[2],p[1]=v[1],p[2]=v[0];    \
    else { if(VAL(v[2])<VAL(v[0])) p[0]=v[1],p[1]=v[2],p[2]=v[0];    \
                              else p[0]=v[1],p[1]=v[0],p[2]=v[2]; }  \
  } else {                                                           \
     if(VAL(v[2])<VAL(v[0]))        p[0]=v[2],p[1]=v[0],p[2]=v[1];   \
     else { if(VAL(v[2])<VAL(v[1])) p[0]=v[0],p[1]=v[2],p[2]=v[1];   \
                               else p[0]=v[0],p[1]=v[1],p[2]=v[2]; } \
  } while(0)
#define MERGE_SORT() \
  do {                                                                 \
    uint i=0, n=An, base=-n, odd=0, c=0, b=1;                          \
    for(;;) {                                                          \
      DATA *p;                                                         \
      if((c&1)==0) {                                                   \
        base+=n, n+=(odd&1), c|=1, b^=1;                               \
        while(n>3) odd<<=1,odd|=(n&1),n>>=1,c<<=1,b^=1;                \
      } else                                                           \
        base-=n-(odd&1),n<<=1,n-=(odd&1),odd>>=1,c>>=1;                \
      if(c==0) break;                                                  \
      p = buf[b]+base;                                                 \
      if(n==2) {                                                       \
        DATA v[2]; SETVAL(v[0],i), SETVAL(v[1],i+1);                   \
        MERGE_2(p,v);                                                  \
        i+=2;                                                          \
      } else if(n==3) {                                                \
        DATA v[3]; SETVAL(v[0],i), SETVAL(v[1],i+1), SETVAL(v[2],i+2); \
        MERGE_3(p,v);                                                  \
        i+=3;                                                          \
      } else {                                                         \
        const uint na = n>>1, nb = (n+1)>>1;                           \
        const DATA *ap = buf[b^1]+base, *ae = ap+na;                   \
        DATA *bp = p+na, *be = bp+nb;                                  \
        for(;;) {                                                      \
          if(VAL((*bp))<VAL((*ap))) {                                  \
            *p++=*bp++;                                                \
            if(bp!=be) continue;                                       \
            do *p++=*ap++; while(ap!=ae);                              \
            break;                                                     \
          } else {                                                     \
            *p++=*ap++;                                                \
            if(ap==ae) break;                                          \
          }                                                            \
        }                                                              \
      }                                                                \
    }                                                                  \
  } while(0)

static void merge_sortv(T *out, const T *A, uint An, const unsigned stride,
                        T *work)
{
  T *buf[2]={out,work};
#define DATA T
#define VAL(x) x
#define SETVAL(x,ai) x=*A,INC_PTR(A,stride)
  MERGE_SORT();
#undef SETVAL
#undef VAL
#undef DATA
}

static void merge_sortp0(uint *idx, const T *A, uint An, const unsigned stride,
                         sort_data *work)
{
  sort_data *buf[2]={work+An,work};
#define DATA sort_data
#define VAL(x) x.v
#define SETVAL(x,ai) x.v=*A,INC_PTR(A,stride),x.i=ai
  MERGE_SORT();
#undef SETVAL
#undef VAL
#undef DATA
  {
    const sort_data *p = buf[0], *pe = p+An;
    do *idx++ = (p++)->i; while(p!=pe);
  }
}

static void merge_sortp(uint *idx, const T *A, uint An, const unsigned stride,
                        sort_data *work)
{
  sort_data *buf[2]={work+An,work};
#define DATA sort_data
#define VAL(x) x.v
#define SETVAL(x,ai) x.i=idx[ai],x.v=INDEX_PTR(A,stride,x.i)
  MERGE_SORT();
#undef SETVAL
#undef VAL
#undef DATA
  {
    const sort_data *p = buf[0], *pe = p+An;
    do *idx++ = (p++)->i; while(p!=pe);
  }
}

#undef MERGE_SORT
#undef MERGE_3
#undef MERGE_2

/*------------------------------------------------------------------------------
  
  Hybrid Stable Sort
  
  low-overhead merge sort when n is small,
  otherwise asymptotically superior radix sort

  result = O(n) sort with good performance for all n
  
  A, n, stride : specifices the input, stride in bytes
  out : the sorted values or indices on output

  ----------------------------------------------------------------------------*/

void sortv(T *out, const T *A, uint n, unsigned stride,
           buffer *buf, int resize)
{
  if(n<2) {
    if(n==0) return;
    *out = *A;
  } else {
    T *w;
    uint work_size = align_as(T,buf->n+n*sizeof(T));
    if(resize) buffer_reserve(buf,work_size);
    else if(buf->max<work_size)
      failwith(__FILE__ ": sortv() not given enough workspace");
    w = align_ptr(T,buf->ptr,buf->n);
    if(n<DIGIT_VALUES) merge_sortv(out, A,n,stride, w);
    else               radix_sortv(out, A,n,stride, w);
  }
}

void sortp(uint *out, int start_perm, const T *A, uint n, unsigned stride,
           buffer *buf, int resize)
{
  if(n<2) {
    if(n==0) return;
    *out = 0;
  } else {
    sort_data *w;
    uint work_size = align_as(sort_data,buf->n+n*2*sizeof(sort_data));
    if(resize) buffer_reserve(buf,work_size);
    else if(buf->max<work_size)
      failwith(__FILE__ ": sortp() not given enough workspace");
    w = align_ptr(sort_data,buf->ptr,buf->n);
    if(n<DIGIT_VALUES)
      if(start_perm) merge_sortp (out, A,n,stride, w);
      else           merge_sortp0(out, A,n,stride, w);
    else
      radix_sortp(out,start_perm, A,n,stride, w);
  }
}

#undef DIGIT_BITS
#undef DIGIT_VALUES
#undef DIGIT_MASK
#undef CEILDIV
#undef DIGITS
#undef VALUE_BITS
#undef COUNT_SIZE

#undef INDEX_PTR
#undef INC_PTR

#undef sortp
#undef sortv

#undef merge_sortp
#undef merge_sortp0
#undef merge_sortv
#undef radix_sortp
#undef radix_passp_be
#undef radix_passp0_be
#undef radix_passp_e
#undef radix_passp_m
#undef radix_passp_b
#undef radix_passp0_b
#undef radix_sortv
#undef radix_passv
#undef radix_zeros
#undef radix_offsets
#undef radix_count
#undef sort_data

