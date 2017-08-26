/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaAsmTest
 Version: 1.0
 Description: Perform various tests on the operaFITSImage class.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaFITSImage.h"

/*! \file operaAsmTest.cpp */

using namespace std;

/*! 
 * operaAsmTest
 * \author Doug Teeple
 * \ingroup test
 */
#if 0
operaFITSImage& operator-=(operaFITSImage* b) {
	 float *p = (float *)pixptr; 
	 float *bp = (float *)b->pixptr; 
	 unsigned n = npixels; 
	 while (n--) *p++ -= *bp++;
	 if (b->istemp) delete b;
	 return *this;
 };
 /*
 * foo is 93 instructions MacOSX...
 */
 opera:test teeple$ otool -Vvt operaAsmTest.o
 operaAsmTest.o:
 (__TEXT,__text) section
 __Z3fooP14operaFITSImageS0_:
 0000000000000000	pushq	%rbp
 0000000000000001	movq	%rsp,%rbp
 0000000000000004	pushq	%rbx
 0000000000000005	subq	$0x08,%rsp
 0000000000000009	movq	%rsi,%rbx
 000000000000000c	movq	0x50(%rdi),%rax
 0000000000000010	movq	0x50(%rsi),%rdx
 0000000000000014	movl	0x40(%rdi),%esi
 0000000000000017	testl	%esi,%esi
 0000000000000019	je	0x0000003a
 000000000000001b	xorl	%ecx,%ecx
 000000000000001d	nopl	(%rax)
 0000000000000020	movss	(%rax),%xmm0			<-- top of assign loop
 0000000000000024	subss	(%rdx),%xmm0			| 2
 0000000000000028	movss	%xmm0,(%rax)			| 3
 000000000000002c	addq	$0x04,%rax				| 4
 0000000000000030	addq	$0x04,%rdx				| 5
 0000000000000034	incl	%ecx					| 6
 0000000000000036	cmpl	%esi,%ecx				| 7
 0000000000000038	jne	0x00000020					<-- bottom of assign loop
 000000000000003a	cmpb	$__Z3fooP14operaFITSImageS0_,0x4c(%rbx)
 000000000000003e	jne	0x00000047
 0000000000000040	addq	$0x08,%rsp
 0000000000000044	popq	%rbx
 0000000000000045	leave
 0000000000000046	ret
 0000000000000047	movq	%rbx,%rdi
 000000000000004a	callq	__ZN14operaFITSImageD1Ev
 000000000000004f	movq	%rbx,%rdi
 0000000000000052	addq	$0x08,%rsp
 0000000000000056	popq	%rbx
 0000000000000057	leave
 0000000000000058	jmp	__ZdlPv
 000000000000005d	nopl	(%rax)
 /*
  * Linux...
  * 6 instructions...
 */
void foo(operaFITSImage *a, operaFITSImage *b) {
	0:	55                   	push   %ebp
	1:	89 e5                	mov    %esp,%ebp
	3:	56                   	push   %esi
	4:	53                   	push   %ebx
	5:	83 ec 10             	sub    $0x10,%esp
	8:	8b 45 08             	mov    0x8(%ebp),%eax
b:	8b 75 0c					mov    0xc(%ebp),%esi

	operaFITSImage& operator-=(operaFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
	e:	8b 50 2c             	mov    0x2c(%eax),%edx
		operaFITSImage& operator-=(operaFITSImage& b) {
			float *p = (float *)pixptr; 
			11:	8b 48 3c             	mov    0x3c(%eax),%ecx
			float *bp = (float *)b.pixptr; 
			14:	8b 5e 3c             	mov    0x3c(%esi),%ebx
			unsigned n = npixels; 
			while (n--) *p++ -= *bp++;
			17:	85 d2                	test   %edx,%edx							
			19:	74 16                	je     31 <_Z3fooP14operaFITSImageS0_+0x31>
			1b:	31 c0                	xor    %eax,%eax
			1d:	8d 76 00             	lea    0x0(%esi),%esi
			20:	d9 04 01             	flds   (%ecx,%eax,1)						<-- top of assign loop
			23:	d8 24 03             	fsubs  (%ebx,%eax,1)						| 2
			26:	d9 1c 01             	fstps  (%ecx,%eax,1)						| 3
			29:	83 c0 04             	add    $0x4,%eax							| 4
			2c:	83 ea 01             	sub    $0x1,%edx							| 5
			2f:	75 ef                	jne    20 <_Z3fooP14operaFITSImageS0_+0x20>	<-- bottom of assign loop
			if (b.istemp) delete &b;
			31:	80 7e 38 00          	cmpb   $0x0,0x38(%esi)
			35:	74 19                	je     50 <_Z3fooP14operaFITSImageS0_+0x50>
			37:	89 34 24             	mov    %esi,(%esp)
			3a:	e8 fc ff ff ff       	call   3b <_Z3fooP14operaFITSImageS0_+0x3b>
			3f:	89 75 08             	mov    %esi,0x8(%ebp)
			*a -= *b;
		}
		42:	83 c4 10             	add    $0x10,%esp
		45:	5b                   	pop    %ebx
		46:	5e                   	pop    %esi
		47:	5d                   	pop    %ebp
		48:	e9 fc ff ff ff       	jmp    49 <_Z3fooP14operaFITSImageS0_+0x49>
		4d:	8d 76 00             	lea    0x0(%esi),%esi
		50:	83 c4 10             	add    $0x10,%esp
		53:	5b                   	pop    %ebx
		54:	5e                   	pop    %esi
		55:	5d                   	pop    %ebp
		56:	c3                   	ret    
		57:	89 f6                	mov    %esi,%esi
		59:	8d bc 27 00 00 00 00 	lea    0x0(%edi,%eiz,1),%edi
		
#endif

void foo(operaFITSImage *a, operaFITSImage *b) {
	*a -= *b;
}
int main()
{
	operaFITSImage a("a.fits", true);
	operaFITSImage b("b.fits", true);
	a -= b;
	
	return 0;
}
