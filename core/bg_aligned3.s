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

.file "bg_aligned3.s"

.globl bg_aligned3
.type  bg_aligned3, @function
.size  bg_aligned3, 48

.section ".text"
.align 2

bg_aligned3:
  andi.    r0,r3,15
  clrlwi   r9,r4,28
  cmpwi    cr7,r9,0
  li       r3,0
  li       r0,0
  bne-     .L.das_label.58
  andi.    r9,r5,15
  bne-     cr7,.L.das_label.58
  bne-     .L.das_label.58
  li       r0,1
 .L.das_label.58:
  stw      r0,0(r6)
  blr     


.ident "GCC: (GNU) 3.2"
