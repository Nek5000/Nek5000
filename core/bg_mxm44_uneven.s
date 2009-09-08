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

.file "bg_mxm44_uneven.s"

.globl bg_mxm44_uneven
.type  bg_mxm44_uneven, @function
.size  bg_mxm44_uneven, 220

.section ".text"
.align 2

bg_mxm44_uneven:
  stwu     r1,-64(r1)
  mflr     r0
  stw      r0,68(r1)
  stw      r28,48(r1)
  lwz      r9,0(r8)
  addi     r8,r1,8
  lwz      r28,0(r4)
  srawi    r0,r9,2
  addze    r0,r0
  stw      r23,28(r1)
  rlwinm   r0,r0,2,0,29
  stw      r24,32(r1)
  subf     r0,r0,r9
  stw      r25,36(r1)
  subf     r9,r0,r9
  stw      r0,12(r1)
  stw      r9,8(r1)
  mr       r23,r4
  stw      r26,40(r1)
  mr       r24,r6
  stw      r27,44(r1)
  mr       r26,r5
  stw      r29,52(r1)
  mr       r27,r7
  mr       r25,r3
  lwz      r29,0(r6)
  crclr    4*cr1+eq
  bl       bg_mxm44
  addi     r8,r1,12
  lwz      r0,8(r1)
  mr       r3,r25
  mr       r4,r23
  mr       r6,r24
  mullw    r28,r28,r0
  mullw    r29,r29,r0
  rlwinm   r28,r28,3,0,28
  add      r27,r27,r28
  mr       r7,r27
  rlwinm   r29,r29,3,0,28
  add      r26,r26,r29
  mr       r5,r26
  crclr    4*cr1+eq
  bl       mxm44_0
  lwz      r29,52(r1)
  lwz      r23,28(r1)
  mr       r3,r0
  lwz      r24,32(r1)
  lwz      r0,68(r1)
  lwz      r25,36(r1)
  lwz      r26,40(r1)
  mtlr     r0
  lwz      r27,44(r1)
  lwz      r28,48(r1)
  addi     r1,r1,64
  blr     

.ident "GCC: (GNU) 3.2"
