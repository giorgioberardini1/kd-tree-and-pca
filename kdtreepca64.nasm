; ---------------------------------------------------------
; PageRank con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software
; installabili mrdiante il packaging tool del sistema
; operativo; per esempio, su Ubuntu, mrdiante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf64 pagerank64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		2.3
;
;align 32
;vec1:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 32
;vec2:		resq	4

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in rax
; l'indirizzo del primo bytes del blocco allocato
; (funziona mrdiante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mrdiante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block
default rel
%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------

input equ 8

indice equ 0

global prova

msg	db 'n:',0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		; rdi = indirizzo della struct input

		; esempio: stampa input->n e di input->k
		; rdi contiente l'indirizzo della struttura contenente i parametri
		; [rdi] contiene l'indirizzo della stringa con il nome del file
		; [rdi+8] contiene l'indirizzo di partenza del data set
		; [rdi+16] contiene l'indirizzo di partenza del query set
		;movsx rax, dword[rdi+24]		; [rdi+16] contiene n
		;prints msg
		;printi rax
		;prints nl
		;movsx rax, dword[rdi+28]		; a 4 byte da n si trova k
		;printi rax
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante



global euclideanDistanceAssembly

euclideanDistanceAssembly:


				push		rbp							; salva il Base Pointer
				mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
				pushaq									; salva i registri generali


				; rdi   ;punto p
				; rsi   ;punto q
				; rdx   ;k*4
				; rcx   ;distance


				mov r8, 0 ; indice i



				vxorps ymm7,ymm7 ;registro risultato
				vxorps ymm6,ymm6 ;registro di supporto a versione scalare


cicloVect:
					cmp rdx,32
					jl cicloResto
					vmovaps ymm0,[rdi+r8]  ;prendo i p[i]..p[i+8]

					vsubps ymm0,[rsi+r8]   ;prendo i q[i]..q[i+8] e li sottraggo ai p
					;vmovss [rcx], xmm0
					;mov rdx, 0
					vmulps ymm0,ymm0       ;quadrato
					vaddps ymm7,ymm0       ;salvo in ymm7 le somme dei quadrati


					add r8,32
					sub rdx,32

					jmp cicloVect


cicloResto:
					cmp rdx,4              ;verifica se hai almeno un elemento
					jl finalize
					vmovss xmm0, [rdi+r8]  ;prendo i p[i]..p[i+4]
			   		vsubss xmm0,[rsi+r8]    ;prendo i q[i]..q[i*4] e li sottraggo ai p
					vmulss xmm0,xmm0         ;quadrato
					vaddss xmm6,xmm0         ;salvo in ymm7 le somme dei quadrati

					add r8,4
					sub rdx,4

					jmp cicloResto



finalize:
					vhaddps ymm7,ymm7 ;sommo con parallel reduction
					vhaddps ymm7,ymm7 ;sommo con parallel reduction
					;vhaddps ymm7,ymm7 ;sommo con parallel reduction

					;vxorps ymm1,ymm1
					vperm2f128 ymm1, ymm7, ymm7, 10000001b
					vaddss  xmm1,xmm7
					vaddss  xmm1,xmm6
					vsqrtss xmm1,xmm1   ;radice quadrata

					vmovss [rcx], xmm1



					popaq						; ripristina i registri generali
					mov		rsp, rbp			; ripristina lo Stack Pointer
					pop		rbp					; ripristina il Base Pointer
					ret							; torna alla funzione C chiamante






global scalarProductAssembly

scalarProductAssembly:

		push		rbp							; salva il Base Pointer
		mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
		pushaq									; salva i registri generali

		; rdi   ;v1
		; rsi   ;n
		; rdx   ;res

		mov r8, 0  ; indice i

		vxorps ymm4,ymm4
		vxorps ymm5,ymm5
		vxorps ymm6,ymm6
		vxorps ymm7,ymm7


cicloUnroll64PS:

		cmp rsi, 128
		jl cicloUnroll16PS

		vmovaps ymm0, [rdi+r8]
		vmovaps ymm1, [rdi+r8+32]
		vmovaps ymm2, [rdi+r8+64]
		vmovaps ymm3, [rdi+r8+96]

		vmulps ymm0, ymm0   ;esegui i quadrati
		vmulps ymm1, ymm1
		vmulps ymm2, ymm2
		vmulps ymm3, ymm3

		vaddps ymm4,ymm0  ;accumula
		vaddps ymm5,ymm1
		vaddps ymm6,ymm2
		vaddps ymm7,ymm3


		add r8, 128
		sub rsi,128
		jmp cicloUnroll64PS

cicloUnroll16PS:

			cmp rsi, 32
			jl finePS
			vmovaps ymm0, [rdi+r8]      ;prendo i v1[i]..v1[i+4]
			vmulps ymm0, ymm0   	    ;moltiplico per se stesso
			vaddps ymm4,ymm0		    ;sommo nel registro xmm4
			add r8, 32
			sub rsi,32
			jmp cicloUnroll16PS

finePS:

		;accumula tutto in xmm7


		vaddps ymm7,ymm4
		vaddps ymm7,ymm5
		vaddps ymm7,ymm6


	    vhaddps ymm7,ymm7 ; sommo
	    vhaddps ymm7,ymm7 ; sommo

		vperm2f128 ymm1, ymm7, ymm7, 10000001b

		vaddss  xmm1,xmm7

		vmovss [rdx],xmm1


		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante



global centerDatasetAssembly

centerDatasetAssembly:

		push		rbp							; salva il Base Pointer
		mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
		pushaq									; salva i registri generali




			; 	rdi <-> eax dataset
			; 	rsi <-> ebx n*4
			;	rdx	k*4
			;	rcx	medieDS

			; r8  <-> ecx n*4

			; xmm0  <->  n


			mov r8,rsi
			vmovss xmm6,xmm0
			vxorps ymm7,ymm7   ;accumulatore vettori
			vxorps ymm9, ymm9
			mov r9,0


cicloMedie64:
			cmp r8, 128
			jl cicloMedie16

			vmovaps ymm0,[rdi+r9]	 ;dataset
			vaddps  ymm0,[rdi+r9+32]
			vaddps  ymm0,[rdi+r9+64]
			vaddps  ymm0,[rdi+r9+96]
			vaddps  ymm7,ymm0

			add r9,128
			sub r8,128
			jmp cicloMedie64



cicloMedie16:
			cmp r8, 32
			jl finaleMedie
			vaddps ymm7, [rdi+r9]

			add r9, 32
			sub r8, 32
			jmp cicloMedie16


finaleMedie:


			vhaddps ymm7,ymm7
			vhaddps ymm7,ymm7

			vperm2f128 ymm1, ymm7, ymm7, 10000001b

			vaddss xmm7,xmm1

			vdivss xmm7,xmm6

			vmovss [rcx],xmm7


			;in xmm7 parte bassa abbiamo la media


			mov r9,0
			vshufps xmm7,xmm7, 00000000
			vperm2f128 ymm9, ymm7, ymm7, 00000000b
			vmovaps ymm7,ymm9

														; vshufps ymm7,ymm7, 00000000
														;vbroadcastf128 ymm7,[rcx]


			; rdi <-> eax dataset
			; rsi <-> ebx n*4
			; r8  <-> ecx n*4
			; rdx <-> edx medieDs
			; xmm0  <->  n

cicloAggDS64:
				cmp rsi, 128
				jl cicloAggDS16

 				vmovaps ymm0,[rdi+r9]
 				vmovaps ymm1,[rdi+r9+32]
 				vmovaps ymm2,[rdi+r9+64]
 				vmovaps ymm3,[rdi+r9+96]

				vsubps ymm0,ymm7
				vsubps ymm1,ymm7
				vsubps ymm2,ymm7
				vsubps ymm3,ymm7

				vmovaps [rdi+r9],ymm0 ;aggiornamento ds
				vmovaps [rdi+r9+32],ymm1 ;aggiornamento ds
				vmovaps [rdi+r9+64],ymm2 ;aggiornamento ds
				vmovaps [rdi+r9+96],ymm3 ;aggiornamento ds

				add r9, 128
				sub rsi, 128

				jmp cicloAggDS64


 cicloAggDS16:  ;si sottrae a ogni elemento (scandito per colonna) la media

				cmp rsi, 32
				jl fineDS

 				vmovaps ymm0,[rdi+r9]
				vsubps ymm0,ymm7
				vmovaps [rdi+r9],ymm0 ;aggiornamento ds

				add r9,32
				sub rsi, 32

				jmp cicloAggDS16

	fineDS:



		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante










global queryPointASMNoPCA




queryPointASMNoPCA:
			push		rbp							; salva il Base Pointer
			mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
			pushaq									; salva i registri generali


			;rdi = eax  u
			;rsi = ebx  ds
			;rdx = edi  k*4
            mov r8, 0

cicloQueryPoint64NoPCA:

			cmp rdx, 128
			jl cicloQueryPoint16NoPCA

			vmovaps ymm0, [rsi+r8]
			vmovaps ymm1, [rsi+r8+32]
			vmovaps ymm2, [rsi+r8+64]
			vmovaps ymm3, [rsi+r8+96]

			vmovaps [rdi+r8], ymm0
			vmovaps [rdi+r8+32], ymm1
			vmovaps [rdi+r8+64], ymm2
			vmovaps [rdi+r8+96], ymm3

			add r8, 64
			sub rdx, 64

			jmp cicloQueryPoint64NoPCA


cicloQueryPoint16NoPCA:

            cmp rdx, 32
            jl fineQPNoPCA

            vmovaps ymm0, [rsi+r8]		;QS
            vmovaps [rdi+r8], ymm0

            add r8, 32
            sub rdx, 32
            jmp cicloQueryPoint16NoPCA


fineQPNoPCA:

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante





global normalizeVAssembly

normalizeVAssembly:
			push		rbp							; salva il Base Pointer
			mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
			pushaq									; salva i registri generali


			;rdi = v
			;rsi = 4*k
			;xmm0 = normV
			mov r8, 0
			vsqrtss xmm0, xmm0
			vshufps ymm0, ymm0, 0 ;in xmm0 ho la norma
			vperm2f128 ymm0, ymm0, ymm0, 0 ;propago la norma su tutto ymm0


cicloUnrollNormalizza128:

			cmp rsi, 128
			jl cicloNormalizza32

			vmovaps ymm1, [rdi+r8]			;accedo al vettore da normalizzare
			vmovaps ymm2, [rdi+r8+32]
			vmovaps ymm3, [rdi+r8+64]
			vmovaps ymm4, [rdi+r8+96]

			vdivps ymm1, ymm0						;divido per la norma
			vdivps ymm2, ymm0
			vdivps ymm3, ymm0
			vdivps ymm4, ymm0

			vmovaps [rdi+r8], ymm1
			vmovaps [rdi+r8+32], ymm2		;riscrivo sul vettore il valore normalizzato
			vmovaps [rdi+r8+64], ymm3
			vmovaps [rdi+r8+96], ymm4

			sub rsi, 128
			add r8, 128

			jmp cicloUnrollNormalizza128

cicloNormalizza32:

			cmp rsi, 32
			jl fineNormalizza

			vmovaps ymm1, [rdi+r8]

			vdivps ymm1, ymm0

			vmovaps [rdi+r8], ymm1

			sub rsi, 32
			add r8, 32

			jmp cicloNormalizza32

fineNormalizza:

			popaq						; ripristina i registri generali
			mov		rsp, rbp			; ripristina lo Stack Pointer
			pop		rbp					; ripristina il Base Pointer
			ret							; torna alla funzione C chiamante








global resetVectorAssembly


resetVectorAssembly:


			push		rbp							; salva il Base Pointer
			mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
			pushaq									; salva i registri generali


			;rdi		u
			;rsi		n*4

			mov r8, 0
			vxorps ymm0, ymm0

cicloresetVectorAssembly128:

			cmp rsi, 128
			jl cicloresetVectorAssembly32

			vmovaps[rdi+r8], ymm0
			vmovaps[rdi+r8+32], ymm0
			vmovaps[rdi+r8+64], ymm0
			vmovaps[rdi+r8+96], ymm0

			add r8, 128
			sub rsi, 128

			jmp cicloresetVectorAssembly128

cicloresetVectorAssembly32:

			cmp rsi, 32
			jl fineresetVectorAssembly

			vmovaps [rdi+r8], ymm0

			add r8, 32
			sub rsi, 32

			jmp cicloresetVectorAssembly32

fineresetVectorAssembly:




			popaq						; ripristina i registri generali
			mov		rsp, rbp			; ripristina lo Stack Pointer
			pop		rbp					; ripristina il Base Pointer
			ret							; torna alla funzione C chiamante


global updateU


updateU:


			push		rbp							; salva il Base Pointer
			mov			rbp, rsp					; il Base Pointer punta al Record di Attivazione corrente
			pushaq									; salva i registri generali



			;rdi = eax 		dataset
			;rsi = ebx 		u
			;xmm0 = xmm6	v[p]
			;rdx =edx		n*4

			;r8= esi


			mov r8, 0
			vxorps ymm9, ymm9

			vshufps xmm0,xmm0, 00000000
			vperm2f128 ymm9, ymm0, ymm0, 00000000b


cicloColonneDataset64:

			cmp rdx,256
			jl cicloColonneDataset16

			vmovaps ymm0, [rdi+r8]  ;leggi 4 elem da dataset
			vmovaps ymm1, [rdi+r8+32]
			vmovaps ymm2, [rdi+r8+64]
			vmovaps ymm3, [rdi+r8+96]
			vmovaps ymm4, [rdi+r8+128]
			vmovaps ymm5, [rdi+r8+160]
			vmovaps ymm6, [rdi+r8+192]
			vmovaps ymm7, [rdi+r8+224]
			

			vmulps ymm0, ymm9   ;effettua prodotto ds*v
			vmulps ymm1, ymm9
			vmulps ymm2, ymm9
			vmulps ymm3, ymm9
			vmulps ymm4, ymm9
			vmulps ymm5, ymm9
			vmulps ymm6, ymm9
			vmulps ymm7, ymm9

			vaddps ymm0, [rsi+r8]   ;aggiungi a u[i]
			vaddps ymm1, [rsi+r8+32]
			vaddps ymm2, [rsi+r8+64]
			vaddps ymm3, [rsi+r8+96]
			vaddps ymm4, [rsi+r8+128]   ;aggiungi a u[i]
			vaddps ymm5, [rsi+r8+160]
			vaddps ymm6, [rsi+r8+192]
			vaddps ymm7, [rsi+r8+224]

			vmovaps [rsi+r8],ymm0   ;scrivi u[i] in memoria
			vmovaps [rsi+r8+32],ymm1   ;scrivi u[i] in memoria
			vmovaps [rsi+r8+64],ymm2   ;scrivi u[i] in memoria
			vmovaps [rsi+r8+96],ymm3   ;scrivi u[i] in memoria
			vmovaps [rsi+r8+128],ymm4
			vmovaps [rsi+r8+160],ymm5
			vmovaps [rsi+r8+192],ymm6
			vmovaps [rsi+r8+224],ymm7

			add r8,256
			sub rdx,256
			jmp cicloColonneDataset64

cicloColonneDataset16:

			cmp rdx,32
			jl fineColonneDataset

			vmovaps ymm0, [rdi+r8]  ;leggi 4 elem da dataset
			vmulps  ymm0,ymm9        ;v[p]*dataset[i..4]
			vaddps ymm0, [rsi+r8]	;leggi 4 elem da u e aggiungi
			vmovaps [rsi+r8],ymm0   ;scrivi u[i] in memoria

			add r8,32
			sub rdx,32
			jmp cicloColonneDataset16

fineColonneDataset:


			popaq						; ripristina i registri generali
			mov		rsp, rbp			; ripristina lo Stack Pointer
			pop		rbp					; ripristina il Base Pointer
			ret							; torna alla funzione C chiamante



global prodottoDsUV2

prodottoDsUV2:


			push		rbp							; salva il Base Pointer
			mov			rbp, rsp					; il Base Pointer punta al Record di Attivazione corrente
			pushaq									; salva i registri generali


			;rdi = eax		dataset
			;rsi = ebx		u
			;rdx = edx		n*4
			;rcx			v
			;xmm0
			;r8 = esi
			;r9 = ecx
			;r10 = edi

			mov r8, 0
			mov r9, 0
			mov r10, rdx
			mov r11, rdx 		;r11=n*4

			vmovss xmm5, xmm0


			vxorps ymm1,ymm1 ;accumlatore v[0]
			vxorps ymm2,ymm2 ;accumlatore v[1]
			vxorps ymm3,ymm3 ;accumlatore v[2]
			vxorps ymm4,ymm4 ;accumlatore v[3]

			
cicloPrimaColonnaProdottoDsUV2256:
			
			cmp r11, 256
			jl cicloPrimaColonnaProdottoDsUV232
			
			vmovaps ymm8, [rdi+r8]		;ds
			vmovaps ymm9, [rdi+r8+32]
			vmovaps ymm10, [rdi+r8+64]
			vmovaps ymm11, [rdi+r8+96]
			vmovaps ymm12, [rdi+r8+128]
			vmovaps ymm13, [rdi+r8+160]
			vmovaps ymm14, [rdi+r8+192]
			vmovaps ymm15, [rdi+r8+224]
			
			vmulps  ymm8, [rsi+r9]		;ds*u
			vmulps  ymm9, [rsi+r9+32]
			vmulps  ymm10, [rsi+r9+64]
			vmulps  ymm11, [rsi+r9+96]
			vmulps  ymm12, [rsi+r9+128]
			vmulps  ymm13, [rsi+r9+160]
			vmulps  ymm14, [rsi+r9+192]
			vmulps  ymm15, [rsi+r9+224]
			
			vaddps  ymm1, ymm8		;v[0]
			vaddps  ymm1, ymm9
			vaddps  ymm1, ymm10
			vaddps  ymm1, ymm11
			vaddps  ymm1, ymm12
			vaddps  ymm1, ymm13
			vaddps  ymm1, ymm14
			vaddps  ymm1, ymm15
			
			add r8, 256
			add r9, 256
			sub r11, 256
			jmp cicloPrimaColonnaProdottoDsUV2256

cicloPrimaColonnaProdottoDsUV232:

			cmp r11, 32
			jl cicloSecondaColonnaProdottoDsUV2Init

			vmovaps ymm0, [rdi+r8]  ;ds[i..4]
			vmulps  ymm0, [rsi+r9]  ;ds[i..4]*u[i..4]
			vaddps  ymm1, ymm0		;v[0]


			add r8,32
			add r9,32
			sub r11, 32
			jmp cicloPrimaColonnaProdottoDsUV232

cicloSecondaColonnaProdottoDsUV2Init:

			;sub r8,32
			add rdx,r10 ;n_padding*2
			mov r9,0
			mov r11, r10 	;n*4
			
cicloSecondaColonnaProdottoDsUV2256:

			cmp r11, 256
			jl cicloSecondaColonnaProdottoDsUV232
			
			vmovaps ymm8, [rdi+r8]		;ds
			vmovaps ymm9, [rdi+r8+32]
			vmovaps ymm10, [rdi+r8+64]
			vmovaps ymm11, [rdi+r8+96]
			vmovaps ymm12, [rdi+r8+128]
			vmovaps ymm13, [rdi+r8+160]
			vmovaps ymm14, [rdi+r8+192]
			vmovaps ymm15, [rdi+r8+224]
			
			vmulps  ymm8, [rsi+r9]		;ds*u
			vmulps  ymm9, [rsi+r9+32]
			vmulps  ymm10, [rsi+r9+64]
			vmulps  ymm11, [rsi+r9+96]
			vmulps  ymm12, [rsi+r9+128]
			vmulps  ymm13, [rsi+r9+160]
			vmulps  ymm14, [rsi+r9+192]
			vmulps  ymm15, [rsi+r9+224]
			
			vaddps  ymm2, ymm8		;v[0]
			vaddps  ymm2, ymm9
			vaddps  ymm2, ymm10
			vaddps  ymm2, ymm11
			vaddps  ymm2, ymm12
			vaddps  ymm2, ymm13
			vaddps  ymm2, ymm14
			vaddps  ymm2, ymm15
			
			add r8, 256
			add r9, 256
			sub r11, 256
			jmp cicloSecondaColonnaProdottoDsUV2256
			

cicloSecondaColonnaProdottoDsUV232:

			cmp r11,32
			jl cicloTerzaColonnaProdottoDsUV2Init

			vmovaps ymm0, [rdi+r8]  ;ds[i..4]
			vmulps  ymm0, [rsi+r9]  ;ds[i..4]*u[i..4]
			vaddps  ymm2, ymm0		;v[1]


			add r8, 32
			add r9, 32
			sub r11, 32
			jmp cicloSecondaColonnaProdottoDsUV232


cicloTerzaColonnaProdottoDsUV2Init:
			
			add rdx,r10 ;n_padding*3 
			mov r9,0
			mov r11, r10 	;n*4
			
cicloTerzaColonnaProdottoDsUV2256:

			cmp r11, 256
			jl cicloTerzaColonnaProdottoDsUV232
			
			vmovaps ymm8, [rdi+r8]		;ds
			vmovaps ymm9, [rdi+r8+32]
			vmovaps ymm10, [rdi+r8+64]
			vmovaps ymm11, [rdi+r8+96]
			vmovaps ymm12, [rdi+r8+128]
			vmovaps ymm13, [rdi+r8+160]
			vmovaps ymm14, [rdi+r8+192]
			vmovaps ymm15, [rdi+r8+224]
			
			vmulps  ymm8, [rsi+r9]		;ds*u
			vmulps  ymm9, [rsi+r9+32]
			vmulps  ymm10, [rsi+r9+64]
			vmulps  ymm11, [rsi+r9+96]
			vmulps  ymm12, [rsi+r9+128]
			vmulps  ymm13, [rsi+r9+160]
			vmulps  ymm14, [rsi+r9+192]
			vmulps  ymm15, [rsi+r9+224]
			
			vaddps  ymm3, ymm8		;v[0]
			vaddps  ymm3, ymm9
			vaddps  ymm3, ymm10
			vaddps  ymm3, ymm11
			vaddps  ymm3, ymm12
			vaddps  ymm3, ymm13
			vaddps  ymm3, ymm14
			vaddps  ymm3, ymm15
			
			add r8, 256
			add r9, 256
			sub r11, 256
			jmp cicloTerzaColonnaProdottoDsUV2256
			

cicloTerzaColonnaProdottoDsUV232:

			cmp r11, 32
			jl cicloQuartaColonnaProdottoDsUV2Init

			vmovaps ymm0, [rdi+r8]  ;ds[i..4]
			vmulps  ymm0, [rsi+r9]  ;ds[i..4]*u[i..4]
			vaddps  ymm3, ymm0		;v[2]


			add r8,32
			add r9,32
			sub r11, 32
			jmp cicloTerzaColonnaProdottoDsUV232


cicloQuartaColonnaProdottoDsUV2Init:

			add rdx,r10 ;n_padding*4
			mov r9,0
			mov r11, r10 	;n*4
		
cicloQuartaColonnaProdottoDsUV2256:

			cmp r11, 256
			jl cicloQuartaColonnaProdottoDsUV232
			
			vmovaps ymm8, [rdi+r8]		;ds
			vmovaps ymm9, [rdi+r8+32]
			vmovaps ymm10, [rdi+r8+64]
			vmovaps ymm11, [rdi+r8+96]
			vmovaps ymm12, [rdi+r8+128]
			vmovaps ymm13, [rdi+r8+160]
			vmovaps ymm14, [rdi+r8+192]
			vmovaps ymm15, [rdi+r8+224]
			
			vmulps  ymm8, [rsi+r9]		;ds*u
			vmulps  ymm9, [rsi+r9+32]
			vmulps  ymm10, [rsi+r9+64]
			vmulps  ymm11, [rsi+r9+96]
			vmulps  ymm12, [rsi+r9+128]
			vmulps  ymm13, [rsi+r9+160]
			vmulps  ymm14, [rsi+r9+192]
			vmulps  ymm15, [rsi+r9+224]
			
			vaddps  ymm4, ymm8		;v[0]
			vaddps  ymm4, ymm9
			vaddps  ymm4, ymm10
			vaddps  ymm4, ymm11
			vaddps  ymm4, ymm12
			vaddps  ymm4, ymm13
			vaddps  ymm4, ymm14
			vaddps  ymm4, ymm15
			
			add r8, 256
			add r9, 256
			sub r11, 256
			jmp cicloQuartaColonnaProdottoDsUV2256


cicloQuartaColonnaProdottoDsUV232:

			cmp r11, 32
			jl fineProdottoDsUV2

			vmovaps ymm0, [rdi+r8]  ;ds[i..4]
			vmulps  ymm0, [rsi+r9]  ;ds[i..4]*u[i..4]
			vaddps  ymm4, ymm0		;v[3]


			add r8, 32
			add r9, 32
			sub r11, 32
			jmp cicloQuartaColonnaProdottoDsUV232

fineProdottoDsUV2:

			vhaddps ymm1,ymm1
			vhaddps ymm2,ymm2
			vhaddps ymm3,ymm3
			vhaddps ymm4,ymm4

			vhaddps ymm1,ymm1
			vhaddps ymm2,ymm2
			vhaddps ymm3,ymm3
			vhaddps ymm4,ymm4

			vmovss xmm0, xmm5  ;prodotto UtU

			vperm2f128 ymm5, ymm1, ymm1, 10000001b
			vaddss xmm1, xmm5
 
			vmovss xmm5, xmm0  	;prodottoUtU sennò lo perdo dopo

			vperm2f128 ymm6, ymm2, ymm2, 10000001b
			vaddss xmm2, xmm6

			vperm2f128 ymm7, ymm3, ymm3, 10000001b
			vaddss xmm3, xmm7

			vperm2f128 ymm0, ymm4, ymm4, 10000001b
			vaddss xmm4, xmm0

			vxorps ymm7,ymm7


			vinsertps xmm7, xmm1, 00000000b
			vinsertps xmm7, xmm2, 00010000b
			vinsertps xmm7, xmm3, 00100000b
			vinsertps xmm7, xmm4, 00110000b


			vshufps xmm5,xmm5, 00000000

			vdivps xmm7,xmm5 ;somma/prodottoUtU


			vmovaps [rcx], xmm7


			popaq						; ripristina i registri generali
			mov		rsp, rbp			; ripristina lo Stack Pointer
			pop		rbp					; ripristina il Base Pointer
			ret							; torna alla funzione C chiamante





global updateMatrix

updateMatrix:

						push		rbp							; salva il Base Pointer
						mov			rbp, rsp					; il Base Pointer punta al Record di Attivazione corrente
						pushaq									; salva i registri generali


						;rdi = eax = U+j*n
						;rsi = ebx = u
						;rdx = ecx = n*4

						mov r8, 0


cicloupdateMatrix128:

						cmp rdx, 128
						jl cicloupdateMatrix32

						vmovaps ymm0, [rsi+r8]
						vmovaps ymm1, [rsi+r8+32]
						vmovaps ymm2, [rsi+r8+64]
						vmovaps ymm3, [rsi+r8+96]

						vmovaps [rdi+r8], ymm0
						vmovaps [rdi+r8+32], ymm1
						vmovaps [rdi+r8+64], ymm2
						vmovaps [rdi+r8+96], ymm3

						add r8, 128
						sub rdx, 128
						jmp cicloupdateMatrix128


cicloupdateMatrix32:

			      cmp rdx, 32
			      jl fineupdateMatrix

			      vmovaps ymm0, [rsi+r8]        ;u[1..4]
			      vmovaps [rdi+r8], ymm0

			      add r8, 32
			      sub rdx, 32
			      jmp cicloupdateMatrix32



fineupdateMatrix:

						popaq						; ripristina i registri generali
						mov		rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante





global updateDataset

updateDataset:

						push		rbp							; salva il Base Pointer
						mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
						pushaq									; salva i registri generali


						;rdi = eax = ds
						;rsi = ebx = u
						;rdx = ecx = v
						;rcx = edx = n*4
						;r8 = edi = k*4
						;r9 = esi = n*k*4


cicloKupdateDatasetVect:

            			mov r10, rcx      ;n*4
						cmp r8,4
						jl fineupdateDataset
						sub r8,4
						vmovss xmm6,[rdx+r8]    	;leggi ultimo elem di v
						vshufps ymm6,ymm6,00000000					
						vperm2f128 ymm6, ymm6, ymm6, 0
						

cicloNupdateDatasetUnrolling64:

						cmp r10,128
						jl cicloNupdateDatasetVect

						sub r10,128     ;n*4 - 8
						sub r9,128      ;n*k - 8

						vmovaps ymm0, [rsi+r10] 	 ;u
						vmovaps ymm1, [rsi+r10+32]
						vmovaps ymm2, [rsi+r10+64]
						vmovaps ymm3, [rsi+r10+96]


						vmulps  ymm0,ymm6	;v*u
						vmulps  ymm1,ymm6
						vmulps  ymm2,ymm6
						vmulps  ymm3,ymm6


						vmovaps ymm4,[rdi+r9]  ;D
						vmovaps ymm5,[rdi+r9+32]
						vmovaps ymm6,[rdi+r9+64]
						vmovaps ymm7,[rdi+r9+96]


						vsubps ymm4,ymm0
						vsubps ymm5,ymm1
						vsubps ymm6,ymm2
						vsubps ymm7,ymm3

						vmovaps [rdi+r9], ymm4 ;scrivi in D
						vmovaps [rdi+r9+32], ymm5 ;scrivi in D
						vmovaps [rdi+r9+64], ymm6 ;scrivi in D
						vmovaps [rdi+r9+96], ymm7

						vmovss xmm6,[rdx+r8]    	;leggi ultimo elem di v
			
						vshufps ymm6,ymm6,00000000
						vperm2f128 ymm6, ymm6, ymm6, 0
						

						jmp cicloNupdateDatasetUnrolling64

cicloNupdateDatasetVect:

						cmp r10,32
						jl cicloKupdateDatasetVect

						sub r10,32
						sub r9,32

						vmovaps ymm0,[rdi+r9]  	;leggi ultimi 4 elem ds
						vmovaps ymm1,[rsi+r10] 		;leggi ultimi 4 elem u
						vmulps ymm1,ymm6	   		;u[1..4]*v[p]
						vsubps ymm0,ymm1 	   		;D-u[1..4]*v[p]

						vmovaps [rdi+r9], ymm0 ;scrivi in D
						jmp cicloNupdateDatasetVect


fineupdateDataset:


						popaq						; ripristina i registri generali
						mov		rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante










global queryPointASM

queryPointASM:

			push		rbp							; salva il Base Pointer
			mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
			pushaq				

		

			;rdi <-> eax   queryPoint
			;rsi <-> ebx   QS 
			;rdx <-> ecx   medieDS 
			;rcx <-> edi   k*4
			;r8  <-> esi		
			
            mov r8, 0		

cicloQueryPoint64:

			cmp rcx, 128
			jl cicloQueryPoint16

			vmovaps ymm0, [rsi+r8]
			vmovaps ymm1, [rsi+r8+32]
			vmovaps ymm2, [rsi+r8+64]
			vmovaps ymm3, [rsi+r8+96]

			vsubps ymm0, [rdx+r8]
			vsubps ymm1, [rdx+r8+32]
			vsubps ymm2, [rdx+r8+64]
			vsubps ymm3, [rdx+r8+96]

			vmovaps [rdi+r8], ymm0
			vmovaps [rdi+r8+32], ymm1
			vmovaps [rdi+r8+64], ymm2
			vmovaps [rdi+r8+96], ymm3

			add r8, 128
			sub rcx, 128
			jmp cicloQueryPoint64


cicloQueryPoint16:

            cmp rcx, 32
            jl fineQP

            vmovaps ymm0, [rsi+r8]		;QS
            vsubps ymm0, [rdx+r8]		; QS-medieDS

            vmovaps [rdi+r8], ymm0

            add r8, 32
            sub rcx, 32
            jmp cicloQueryPoint16



fineQP:


			popaq						; ripristina i registri generali
			mov		rsp, rbp			; ripristina lo Stack Pointer
			pop		rbp					; ripristina il Base Pointer
			ret							; torna alla funzione C chiaman










global productQPV

productQPV:


		
			push		rbp							; salva il Base Pointer
			mov			rbp, rsp				; il Base Pointer punta al Record di Attivazione corrente
			pushaq				

		
		
			;rdi <-> eax   queryPoint
			;rsi <-> ebx   v+i*input->k
			;rdx <-> edi   k*4
			;rcx <-> ecx   queryPointTmp+i
			;r8 <-> esi

            mov r8,0
            vxorps ymm7, ymm7


cicloproductQPV64:

			cmp rdx, 256
			jl cicloproductQPV16

			vmovaps ymm0, [rdi+r8]
			vmovaps ymm1, [rdi+r8+32]
			vmovaps ymm2, [rdi+r8+64]
			vmovaps ymm3, [rdi+r8+96]
			vmovaps ymm4, [rdi+r8+128]
			vmovaps ymm5, [rdi+r8+160]
			vmovaps ymm6, [rdi+r8+192]
			vmovaps ymm8, [rdi+r8+224]

			vmulps ymm0, [rsi+r8]
			vmulps ymm1, [rsi+r8+32]
			vmulps ymm2, [rsi+r8+64]
			vmulps ymm3, [rsi+r8+96]
			vmulps ymm4, [rsi+r8+128]
			vmulps ymm5, [rsi+r8+160]
			vmulps ymm6, [rsi+r8+192]
			vmulps ymm8, [rsi+r8+224]

			vaddps ymm7, ymm0
			vaddps ymm7, ymm1
			vaddps ymm7, ymm2
			vaddps ymm7, ymm3
			vaddps ymm7, ymm4
			vaddps ymm7, ymm5
			vaddps ymm7, ymm6
			vaddps ymm7, ymm8

			add r8, 256
			sub rdx, 256
			jmp cicloproductQPV64

cicloproductQPV16:

           cmp rdx, 32
           jl fineproductQPV

           vmovaps ymm0, [rdi+r8]
           vmulps ymm0, [rsi+r8]
           vaddps ymm7, ymm0

           sub rdx, 32
           add r8, 32

           jmp cicloproductQPV16




fineproductQPV:

			
			vhaddps ymm7, ymm7
			vhaddps ymm7, ymm7	

			vperm2f128 ymm1, ymm7, ymm7, 10000001b

			vaddps ymm7,ymm1

			movss [rcx],xmm7

			popaq						; ripristina i registri generali
			mov		rsp, rbp			; ripristina lo Stack Pointer
			pop		rbp					; ripristina il Base Pointer
			ret							; torna alla funzione C chiaman

           
			