.set r0,0; .set r1,1; .set r2,2; .set r3,3; .set r4,4
.set r5,5; .set r6,6; .set r7,7; .set r8,8; .set r9,9
.set r10,10; .set r11,11; .set r12,12; .set r13,13; .set r14,14
.set r15,15; .set r16,16; .set r17,17; .set r18,18; .set r19,19
.set r20,20; .set r21,21; .set r22,22; .set r23,23; .set r24,24
.set r25,25; .set r26,26; .set r27,27; .set r28,28; .set r29,29
.set r30,30; .set r31,31
.set f0,0; .set f1,1; .set f2,2; .set f3,3; .set f4,4
.set f5,5; .set f6,6; .set f7,7; .set f8,8; .set f9,9
.set f10,10; .set f11,11; .set f12,12; .set f13,13; .set f14,14
.set f15,15; .set f16,16; .set f17,17; .set f18,18; .set f19,19
.set f20,20; .set f21,21; .set f22,22; .set f23,23; .set f24,24
.set f25,25; .set f26,26; .set f27,27; .set f28,28; .set f29,29
.set f30,30; .set f31,31

.file "bg_mxm3.s"

.globl bg_mxm3
.type  bg_mxm3, @function
.size  bg_mxm3, 1412

.section ".text"
.align 2

bg_mxm3:
  stwu     r1,-96(r1)
  mflr     r0
  stw      r0,100(r1)
  andi.    r0,r7,15
  stw      r15,28(r1)
  mr       r15,r8
  stw      r16,32(r1)
  mr       r16,r6
  stw      r25,68(r1)
  mr       r25,r5
  stw      r28,80(r1)
  mr       r28,r4
  stw      r30,88(r1)
  mr       r30,r7
  stw      r31,92(r1)
  mr       r31,r3
  stw      r14,24(r1)
  stw      r17,36(r1)
  stw      r18,40(r1)
  stw      r19,44(r1)
  stw      r20,48(r1)
  stw      r21,52(r1)
  stw      r22,56(r1)
  stw      r23,60(r1)
  stw      r26,72(r1)
  stw      r27,76(r1)
  stw      r29,84(r1)
  bne-     .L.das_label.15
 .L.das_label.3:
  lis      r23,dummy@ha
  li       r12,0
  addi     r23,r23,dummy@l
  li       r13,0
  li       r27,16
  addi     r23,r23,-16
  stfpdux  f13,r23,r27
  stfpdux  f14,r23,r27
  stfpdux  f15,r23,r27
  stfpdux  f16,r23,r27
  stfpdux  f17,r23,r27
  stfpdux  f18,r23,r27
  stfpdux  f19,r23,r27
  stfpdux  f20,r23,r27
  stfpdux  f21,r23,r27
  stfpdux  f22,r23,r27
  stfpdux  f23,r23,r27
  stfpdux  f24,r23,r27
  stfpdux  f25,r23,r27
  stfpdux  f26,r23,r27
  stfpdux  f27,r23,r27
  stfpdux  f28,r23,r27
  stfpdux  f29,r23,r27
  stfpdux  f30,r23,r27
  stfpdux  f31,r23,r27
  lis      r9,10922
  lwz      r4,0(r28)
  ori      r9,r9,43691
  lwz      r15,0(r15)
  lwz      r6,0(r16)
  mulhw    r14,r4,r9
  srawi    r11,r15,31
  srawi    r16,r6,1
  addze    r16,r16
  srawi    r0,r4,31
  mulhw    r15,r15,r9
  subf     r14,r0,r14
  rlwinm   r26,r16,4,0,27
  mulli    r22,r16,80
  subf     r15,r11,r15
  subfic   r22,r22,16
  mfctr    r29
  li       r18,0
  rlwinm   r28,r4,3,0,28
  cmpw     r18,r15
  addi     r28,r28,-32
  addi     r19,r16,-2
  cmpwi    cr7,r12,0
  bge-     .L.das_label.10
 .L.das_label.4:
  mullw    r23,r4,r18
  li       r17,0
  cmpw     r17,r14
  mulli    r23,r23,48
  add      r23,r30,r23
  addi     r23,r23,-16
  bge-     .L.das_label.9
  mullw    r9,r16,r18
  cmpwi    cr6,r19,0
  mulli    r9,r9,96
  mulli    r0,r6,40
  add      r9,r9,r0
  add      r9,r9,r25
  addi     r11,r9,-16
 .L.das_label.5:
  mulli    r20,r17,48
  rlwinm   r0,r4,3,0,28
  add      r20,r31,r20
  subf     r20,r0,r20
  addi     r20,r20,32
  lfpdux   f1,r20,r28
  mr       r21,r11
  lfpdux   f7,r21,r22
  lfpdux   f2,r20,r27
  lfpdux   f3,r20,r27
  lfpdux   f8,r21,r26
  lfpdux   f9,r21,r26
  lfpdux   f4,r20,r28
  lfpdux   f10,r21,r26
  lfpdux   f5,r20,r27
  lfpdux   f11,r21,r26
  lfpdux   f12,r21,r26
  fxpmul   f13,f7,f1
  fxpmul   f16,f8,f1
  lfpdux   f6,r20,r27
  fxpmul   f19,f9,f1
  fxpmul   f22,f10,f1
  fxpmul   f25,f11,f1
  fxpmul   f28,f12,f1
  fxpmul   f14,f7,f2
  lfpdux   f1,r20,r28
  fxpmul   f17,f8,f2
  fxpmul   f20,f9,f2
  fxpmul   f23,f10,f2
  fxpmul   f26,f11,f2
  fxpmul   f29,f12,f2
  fxpmul   f15,f7,f3
  lfpdux   f2,r20,r27
  fxpmul   f18,f8,f3
  fxpmul   f21,f9,f3
  fxpmul   f24,f10,f3
  fxpmul   f27,f11,f3
  fxpmul   f30,f12,f3
  fxcsmadd f13,f7,f4,f13
  lfpdux   f3,r20,r27
  fxcsmadd f16,f8,f4,f16
  fxcsmadd f14,f7,f5,f14
  fxcsmadd f17,f8,f5,f17
  fxcsmadd f19,f9,f4,f19
  fxcsmadd f22,f10,f4,f22
  fxcsmadd f25,f11,f4,f25
  fxcsmadd f28,f12,f4,f28
  fxcsmadd f20,f9,f5,f20
  fxcsmadd f15,f7,f6,f15
  lfpdux   f4,r20,r28
  fxcsmadd f18,f8,f6,f18
  lfpdux   f7,r21,r22
  fxcsmadd f23,f10,f5,f23
  fxcsmadd f26,f11,f5,f26
  fxcsmadd f29,f12,f5,f29
  lfpdux   f8,r21,r26
  fxcsmadd f21,f9,f6,f21
  fxcsmadd f24,f10,f6,f24
  lfpdux   f5,r20,r27
  lfpdux   f9,r21,r26
  lfpdux   f10,r21,r26
  fxcsmadd f27,f11,f6,f27
  fxcsmadd f30,f12,f6,f30
  lfpdux   f11,r21,r26
  beq-     cr6,.L.das_label.6
  mtctr    r19

.__loopk:
  lfpdux   f12,r21,r26
  fxcpmadd f13,f7,f1,f13
  fxcpmadd f16,f8,f1,f16
  lfpdux   f6,r20,r27
  fxcpmadd f19,f9,f1,f19
  fxcpmadd f22,f10,f1,f22
  fxcpmadd f25,f11,f1,f25
  fxcpmadd f28,f12,f1,f28
  fxcpmadd f14,f7,f2,f14
  lfpdux   f1,r20,r28
  fxcpmadd f17,f8,f2,f17
  fxcpmadd f20,f9,f2,f20
  fxcpmadd f23,f10,f2,f23
  fxcpmadd f26,f11,f2,f26
  fxcpmadd f29,f12,f2,f29
  fxcpmadd f15,f7,f3,f15
  lfpdux   f2,r20,r27
  fxcpmadd f18,f8,f3,f18
  fxcpmadd f21,f9,f3,f21
  fxcpmadd f24,f10,f3,f24
  fxcpmadd f27,f11,f3,f27
  fxcpmadd f30,f12,f3,f30
  fxcsmadd f13,f7,f4,f13
  lfpdux   f3,r20,r27
  fxcsmadd f16,f8,f4,f16
  fxcsmadd f14,f7,f5,f14
  fxcsmadd f17,f8,f5,f17
  fxcsmadd f19,f9,f4,f19
  fxcsmadd f22,f10,f4,f22
  fxcsmadd f25,f11,f4,f25
  fxcsmadd f28,f12,f4,f28
  fxcsmadd f20,f9,f5,f20
  fxcsmadd f15,f7,f6,f15
  lfpdux   f4,r20,r28
  fxcsmadd f18,f8,f6,f18
  lfpdux   f7,r21,r22
  fxcsmadd f23,f10,f5,f23
  fxcsmadd f26,f11,f5,f26
  fxcsmadd f29,f12,f5,f29
  lfpdux   f8,r21,r26
  fxcsmadd f21,f9,f6,f21
  fxcsmadd f24,f10,f6,f24
  lfpdux   f5,r20,r27
  lfpdux   f9,r21,r26
  lfpdux   f10,r21,r26
  fxcsmadd f27,f11,f6,f27
  fxcsmadd f30,f12,f6,f30
  lfpdux   f11,r21,r26
  bdnz+    .__loopk
 .L.das_label.6:
  lfpdux   f12,r21,r26
  fxcpmadd f13,f7,f1,f13
  fxcpmadd f14,f7,f2,f14
  fxcpmadd f15,f7,f3,f15
  lfpdux   f6,r20,r27
  fxcpmadd f16,f8,f1,f16
  fxcpmadd f17,f8,f2,f17
  fxcpmadd f18,f8,f3,f18
  fxcsmadd f13,f7,f4,f13
  fxcsmadd f14,f7,f5,f14
  fxcsmadd f15,f7,f6,f15
  fxcsmadd f16,f8,f4,f16
  fxcsmadd f17,f8,f5,f17
  fxcsmadd f18,f8,f6,f18
  stfpdux  f13,r23,r27
  fxcpmadd f19,f9,f1,f19
  fxcpmadd f20,f9,f2,f20
  stfpdux  f14,r23,r27
  fxcpmadd f21,f9,f3,f21
  fxcpmadd f22,f10,f1,f22
  stfpdux  f15,r23,r27
  fxcpmadd f23,f10,f2,f23
  fxcpmadd f24,f10,f3,f24
  stfpdux  f16,r23,r28
  fxcsmadd f19,f9,f4,f19
  fxcsmadd f20,f9,f5,f20
  stfpdux  f17,r23,r27
  fxcsmadd f21,f9,f6,f21
  fxcsmadd f22,f10,f4,f22
  stfpdux  f18,r23,r27
  fxcsmadd f23,f10,f5,f23
  fxcsmadd f24,f10,f6,f24
  stfpdux  f19,r23,r28
  fxcpmadd f25,f11,f1,f25
  fxcpmadd f26,f11,f2,f26
  stfpdux  f20,r23,r27
  fxcpmadd f27,f11,f3,f27
  fxcpmadd f28,f12,f1,f28
  stfpdux  f21,r23,r27
  fxcpmadd f29,f12,f2,f29
  fxcpmadd f30,f12,f3,f30
  stfpdux  f22,r23,r28
  clrlwi   r0,r23,28
  subfic   r10,r13,0
  adde     r9,r10,r13
  neg      r0,r0
  rlwinm   r0,r0,1,31,31
  and.     r10,r0,r9
  beq-     .L.das_label.7
  mr       r13,r23
 .L.das_label.7:
  fxcsmadd f25,f11,f4,f25
  fxcsmadd f26,f11,f5,f26
  stfpdux  f23,r23,r27
  clrlwi   r0,r23,28
  mfcr     r9
  rlwinm   r9,r9,31,31,31
  neg      r0,r0
  rlwinm   r0,r0,1,31,31
  and.     r10,r0,r9
  beq-     .L.das_label.8
  mr       r12,r23
  cmpwi    cr7,r23,0
 .L.das_label.8:
  fxcsmadd f27,f11,f6,f27
  fxcsmadd f28,f12,f4,f28
  stfpdux  f24,r23,r27
  fxcsmadd f29,f12,f5,f29
  fxcsmadd f30,f12,f6,f30
  stfpdux  f25,r23,r28
  stfpdux  f26,r23,r27
  stfpdux  f27,r23,r27
  stfpdux  f28,r23,r28
  stfpdux  f29,r23,r27
  stfpdux  f30,r23,r27
  addi     r17,r17,1
  mulli    r0,r4,40
  cmpw     r17,r14
  subf     r23,r0,r23
  blt+     .L.das_label.5
 .L.das_label.9:
  addi     r18,r18,1
  cmpw     r18,r15
  blt+     .L.das_label.4
 .L.das_label.10:
  mtctr    r29
  lis      r23,dummy@ha
  addi     r23,r23,dummy@l
  addi     r23,r23,-16
  lfpdux   f13,r23,r27
  lfpdux   f14,r23,r27
  lfpdux   f15,r23,r27
  lfpdux   f16,r23,r27
  lfpdux   f17,r23,r27
  lfpdux   f18,r23,r27
  lfpdux   f19,r23,r27
  lfpdux   f20,r23,r27
  lfpdux   f21,r23,r27
  lfpdux   f22,r23,r27
  lfpdux   f23,r23,r27
  lfpdux   f24,r23,r27
  lfpdux   f25,r23,r27
  lfpdux   f26,r23,r27
  lfpdux   f27,r23,r27
  lfpdux   f28,r23,r27
  lfpdux   f29,r23,r27
  lfpdux   f30,r23,r27
  lfpdux   f31,r23,r27
  bne-     cr7,.L.das_label.14
 .L.das_label.11:
  cmpwi    r13,0
  bne-     .L.das_label.13
 .L.das_label.12:
  lwz      r0,100(r1)
  li       r3,3
  lwz      r14,24(r1)
  lwz      r15,28(r1)
  mtlr     r0
  lwz      r16,32(r1)
  lwz      r17,36(r1)
  lwz      r18,40(r1)
  lwz      r19,44(r1)
  lwz      r20,48(r1)
  lwz      r21,52(r1)
  lwz      r22,56(r1)
  lwz      r23,60(r1)
  lwz      r25,68(r1)
  lwz      r26,72(r1)
  lwz      r27,76(r1)
  lwz      r28,80(r1)
  lwz      r29,84(r1)
  lwz      r30,88(r1)
  lwz      r31,92(r1)
  addi     r1,r1,96
  blr     
 .L.das_label.13:
  lis      r3,.rodata@ha
  mr       r4,r13
  addi     r3,r3,.rodata@l
  crclr    4*cr1+eq
  bl       printf
  b        .L.das_label.12
 .L.das_label.14:
  lis      r3,.rodata+0x00000018@ha
  mr       r4,r12
  addi     r3,r3,.rodata+0x00000018@l
  crclr    4*cr1+eq
  bl       printf
  b        .L.das_label.11
 .L.das_label.15:
  lis      r3,.rodata+0x00000030@ha
  mr       r4,r7
  addi     r3,r3,.rodata+0x00000030@l
  crclr    4*cr1+eq
  bl       printf
  b        .L.das_label.3

.comm dummy,512,16


.section ".rodata","a"
.align 2
  .ascii "ALIGNMENT PROBS1: %x\n"
.align 3
  .ascii "ALIGNMENT PROBS: %x\n"
.align 3
  .ascii "C PROB: %x\n"

.ident "GCC: (GNU) 3.2"
