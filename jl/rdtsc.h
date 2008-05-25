
static __inline__ uint64_t rdtsc()
{
  uint64_t x;
  __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
  return x;
}

