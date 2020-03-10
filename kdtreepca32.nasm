; ---------------------------------------------------------
; PQNN con istruzioni SSE a 32 bit
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
; installabili mediante il packaging tool del sistema
; operativo; per esempio, su Ubuntu, mediante i comandi:
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
;     nasm -f elf32 pqnn32.nasm
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
;vec2:		resq	4

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global prova



input		equ		8

x db "ciao",0;


prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione

		; esempio: stampa input->n e di input->k
		mov EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
		movss XMM0, [EAX+12] ;xmm0=8000
		getmem 4,1 ;alloca 4 byte
		mov ECX,1
		add ECX,2
		; [EAX] contiene l'indirizzo della stringa con il nome del file
		; [EAX+4] contiene l'indirizzo di partenza del data set
		; [EAX+8] contiene l'indirizzo di partenza del query set

		;xorps XMM1,XMM1
		;mov [EAX],ECX
		;movss XMM1,[EAX]
		;addss XMM1,XMM0
		;movss [EAX],XMM1
		;mov EAX,[EAX]











		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante



global euclideanDistanceAssembly

euclideanDistanceAssembly:


		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi


		mov eax,[EBP+input]   ;punto p
		mov ecx,[EBP+input+4] ;punto q
		mov edi,[EBP+input+8] ;k*4

		mov esi, 0  ; indice i


		xorps xmm7,xmm7 ;registro risultato
		xorps xmm6,xmm6 ;registro di supporto


cicloUnroll64:
				cmp edi,64
				jl cicloVect

				movaps xmm0,[eax+esi]     ;16 byte
				movaps xmm1,[eax+esi+16]  ;16 byte
				movaps xmm2,[eax+esi+32]  ;16 byte
				movaps xmm3,[eax+esi+48]  ;16 byte

				subps  xmm0,[ecx+esi]
				subps  xmm1,[ecx+esi+16]
				subps  xmm2,[ecx+esi+32]
				subps  xmm3,[ecx+esi+48]


				mulps xmm0,xmm0
				mulps xmm1,xmm1
				mulps xmm2,xmm2
				mulps xmm3,xmm3

				addps xmm7,xmm0
				addps xmm7,xmm1
				addps xmm7,xmm2
				addps xmm7,xmm3

				add esi,64
				sub edi,64


				jmp cicloUnroll64


cicloVect:  cmp edi,16
			jl cicloResto
			movaps xmm0, [eax+esi]  ;prendo i p[i]..p[i+4]
	        subps xmm0,[ecx+esi]    ;prendo i q[i]..q[i*4] e li sottraggo ai p
			mulps xmm0,xmm0         ;quadrato
			addps xmm7,xmm0         ;salvo in xmm7 le somme dei quadrati
			add esi,16
			sub edi,16

			jmp cicloVect


cicloResto: 
			; cmp edi,4              ;verifica se hai almeno un elemento
			; jl finalize
			; movss xmm0, [eax+esi]  ;prendo i p[i]..p[i+4]
	        ; subss xmm0,[ecx+esi]    ;prendo i q[i]..q[i*4] e li sottraggo ai p
			; mulss xmm0,xmm0         ;quadrato
			; addss xmm6,xmm0         ;salvo in xmm7 le somme dei quadrati

			; add esi,4
			; sub edi,4

			; jmp cicloResto



finalize:
			haddps xmm7,xmm7 ; sommo con parallel reduction
			haddps xmm7,xmm7 ; sommo con parallel reduction

			addps  xmm7,xmm6
			sqrtss xmm7,xmm7 ; radice quadrata

			mov esi, [ebp+input+12]
			movss [esi], xmm7



		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante




;riscritto da Giorgio per evitare conflitti nei registri, l'altra versione può essere eliminata
global scalarProductAssembly

scalarProductAssembly:

		push		ebp
		mov			ebp, esp
		push		ebx
		push		esi
		push		edi


		mov eax,[ebp+input]   ;v1
		mov edx,[ebp+input+4] ;n (lunghezza v1)
		mov ecx,[ebp+input+8]

		mov esi, 0  ; indice i

		xorps xmm4,xmm4
		xorps xmm5,xmm5
		xorps xmm6,xmm6
		xorps xmm7,xmm7


cicloUnroll64PS:

		cmp edx, 64
		jl cicloUnroll16PS

		movaps xmm0, [eax+esi]
		movaps xmm1, [eax+esi+16]
		movaps xmm2, [eax+esi+32]
		movaps xmm3, [eax+esi+48]

		mulps xmm0, xmm0   ;esegui i quadrati
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3

		addps xmm4,xmm0  ;accumula
		addps xmm5,xmm1
		addps xmm6,xmm2
		addps xmm7,xmm3


		add esi, 64
		sub edx, 64
		jmp cicloUnroll64PS

cicloUnroll16PS:

			cmp edx, 16
			jl cicloScalarePSInit
			movaps xmm0, [eax+esi]  ;prendo i v1[i]..v1[i+4]
			mulps xmm0, xmm0   	    ;moltiplico per se stesso
			addps xmm4,xmm0		    ;sommo nel registro xmm4
			add esi,16
			sub edx, 16
			jmp cicloUnroll16PS


cicloScalarePSInit:

		xorps xmm1,xmm1 ;riutilizza xmm1


cicloScalarePS:
		; cmp edx, 4
		; jl finePS
		; movss xmm0, [eax + esi]
		; mulss xmm0, xmm0
		; addss xmm1, xmm0
		; add esi, 4
		; sub edx, 4
		; jmp cicloScalarePS


finePS:

		;accumula tutto in xmm7

		addps xmm7,xmm1   ;nel caso in cui NON entra mai in cicloScalarePS aggiunge 0
		addps xmm7,xmm4
		addps xmm7,xmm5
		addps xmm7,xmm6


	    haddps xmm7,xmm7 ; sommo
	    haddps xmm7,xmm7 ; sommo







		movss [ecx],xmm7


		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante


;versione obsoleta
; global scalarProductAssembly

; scalarProductAssembly:

; 		push		ebp							; salva il Base Pointer
; 		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
; 		push		ebx							; salva i registri da preservare
; 		push		esi
; 		push		edi




; 		mov EAX,[EBP+input]   ;v1
; 		mov EDX,[EBP+input+4] ;n (lunghezza v1)


; 		mov esi, 0  ; indice i

; 		xorps xmm1,xmm1  ;accumulatore
; 		xorps xmm2, xmm2 ;accumulatore per somme scalari
; 		xorps xmm6,xmm6
; 		xorps xmm7,xmm7


; cicloUnroll64PS:
; 		cmp edx, 64
; 		;jl cicloUnroll32PS
; 		jl cicloUnroll16PS

; 		movaps xmm0, [eax+esi]
; 		movaps xmm3, [eax+esi+16]
; 		movaps xmm4, [eax+esi+32]
; 		movaps xmm5, [eax+esi+48]

; 		mulps xmm0, xmm0   ;esegui i quadrati
; 		mulps xmm3, xmm3
; 		mulps xmm4, xmm4
; 		mulps xmm5, xmm5

; 		; addps xmm1, xmm0	;accumula in xmm1
; 		; addps xmm1, xmm3
; 		; addps xmm1, xmm4
; 		; addps xmm1, xmm5

; ; i registri liberi sono xmm6, xmm7
; 		addps xmm6, xmm0
; 		addps xmm7,xmm3
; 		addps xmm1, xmm4
; 		addps xmm6, xmm5


; 		add esi, 64
; 		sub edx, 64
; 		jmp cicloUnroll64PS

; ;non serve
; ; cicloUnroll32PS:
; ; 		cmp edx, 32
; ; 		jl cicloUnroll16PS

; ; 		movaps xmm0, [eax+esi]
; ; 		movaps xmm3, [eax+esi+16]

; ; 		mulps xmm0, xmm0
; ; 		addps xmm1, xmm0


; ; 		mulps xmm3, xmm3
; ; 		addps xmm1, xmm3

; ; 		add esi, 32
; ; 		sub edx, 32
; ; 		jmp cicloUnroll32PS

; cicloUnroll16PS:

; 			cmp edx, 16
; 			jl cicloScalarePS
; 			movaps XMM0, [EAX+ESI]  ;prendo i v1[i]..v1[i+4]
; 			mulps XMM0, XMM0   	    ;moltiplico per se stesso
; 			addps XMM1,XMM0		    ;sommo nel registro XMM1
; 			add ESI,16
; 			sub edx, 16
; 			jmp cicloUnroll16PS

; cicloScalarePS:
; 		cmp edx, 4
; 		jl fine
; 		movss xmm0, [eax + esi]
; 		mulss xmm0, xmm0
; 		addss xmm2, xmm0
; 		add esi, 4
; 		sub edx, 4
; 		jmp cicloScalarePS


; fine:  ;alternativa
; 		addps xmm1, xmm6
; 		addps xmm1, xmm7

; 		addps xmm1, xmm2
; 	    haddps xmm1,xmm1 ; sommo
; 	    haddps xmm1,xmm1 ; sommo



; 		getmem 4,1



; 		movss [eax],xmm1


; 		pop	edi									; ripristina i registri da preservare
; 		pop	esi
; 		pop	ebx
; 		mov	esp, ebp							; ripristina lo Stack Pointer
; 		pop	ebp									; ripristina il Base Pointer
; 		ret										; torna alla funzione C chiamante



global funzioneProvaAssembly


funzioneProvaAssembly:

	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi


   mov eax,[EBP+input]   ;v1
   mov edx,[EBP+input+4] ;n
   xorps xmm0,xmm0
   mov edi,[ebp+input+8]
   movss [edi],xmm0


;    movaps xmm7,[eax]
;    ;shufps xmm7,xmm7, 00000000
; ;    getmem 4,1
; ;    mov [eax],input
; ;    movss xmm0,input; xmm0=1
;    ;mov eax, [ebp+input]

; 	xorps xmm0,xmm0

;    vinsertps xmm7,xmm0, 00100000b
;    mov ecx,[ebp+input]
;    movaps [ecx],xmm7







;    movaps [eax],xmm7






	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante




global centerDatasetAssembly

centerDatasetAssembly:



	push		ebp
	mov			ebp, esp
	push		ebx
	push		esi
	push		edi


			mov eax,[ebp+input]    ;dataset
			mov ebx,[ebp+input+4]  ;n*4
			mov ecx,[ebp+input+4]  ;n*4
			mov edx,[ebp+input+12] ;medieDs
			movss xmm6,[ebp+input+16] ;n


			xorps xmm7,xmm7   ;accumulatore vettori
			xorps xmm2, xmm2  ;accumulatore scalari

			xorps xmm0,xmm0
			xorps xmm1,xmm1
			xorps xmm3,xmm3
			xorps xmm4,xmm4
			mov esi,0


cicloMedie64:
			cmp ecx, 64
			jl cicloMedie16
			movaps xmm0,[eax+esi]	 ;dataset

	    	addps xmm0,[eax+esi+16]
	    	addps xmm0,[eax+esi+32]
	    	addps xmm0,[eax+esi+48]


			addps xmm7,xmm0
			add esi,64

			sub ecx, 64
			jmp cicloMedie64



; cicloMedie32:
; 			cmp ecx, 32
; 			jl cicloMedie16
; 			movaps xmm0,[eax+esi]
; 	    	addps xmm0,[eax+esi+16]
; 			addps xmm7,xmm0
; 			add esi,32

; 			sub ecx, 32



cicloMedie16:
			cmp ecx, 16
			jl cicloScalareMedie
			addps xmm7, [eax+esi]

			add esi, 16
			sub ecx, 16
			jmp cicloMedie16




cicloScalareMedie:
			; cmp ecx, 4
			; jl finaleMedie

			; addss xmm2, [eax+esi]
			; add esi, 4
			; sub ecx, 4
			; jmp cicloScalareMedie


finaleMedie:

			; addps xmm0,xmm1
			; addps xmm4,xmm3

			;addps xmm7,xmm1
			;addps xmm7,xmm3

			;addps xmm7,xmm4
			;addps xmm7,xmm0

			;addps xmm7, xmm2

			haddps xmm7,xmm7
			haddps xmm7,xmm7

			divss xmm7,xmm6

			movss [edx],xmm7


;in xmm7 parte bassa abbiamo la media


			mov esi,0
			shufps xmm7,xmm7, 00000000

cicloAggDS64:
				cmp ebx, 64
				jl cicloAggDS16

 				movaps xmm0,[eax+esi]
 				movaps xmm1,[eax+esi+16]
 				movaps xmm2,[eax+esi+32]
 				movaps xmm3,[eax+esi+48]

				subps xmm0,xmm7
				subps xmm1,xmm7
				subps xmm2,xmm7
				subps xmm3,xmm7

				movaps [eax+esi],xmm0 ;aggiornamento ds

				movaps [eax+esi+16],xmm1 ;aggiornamento ds

				movaps [eax+esi+32],xmm2 ;aggiornamento ds

				movaps [eax+esi+48],xmm3 ;aggiornamento ds

				add esi, 64
				sub ebx, 64

				jmp cicloAggDS64


 cicloAggDS16:  ;si sottrae a ogni elemento (scandito per colonna) la media

				cmp ebx, 16
				jl cicloScalareAggDS

 				movaps xmm0,[eax+esi]
				subps xmm0,xmm7

				movaps [eax+esi],xmm0 ;aggiornamento ds
				add esi,16
				sub ebx, 16

				jmp cicloAggDS16

cicloScalareAggDS:

			cmp ebx, 4
			jl fineDS

			movss xmm0, [eax+esi]
			subss xmm0, xmm7

			movss [eax+esi], xmm0
			add esi, 4
			sub ebx, 4

			jmp cicloScalareAggDS

fineDS:

	pop	edi
	pop	esi
	pop	ebx
	mov	esp, ebp
	pop	ebp
	ret



;centra Dataset fatto da Giorgio, le performance peggiorano con questo unrolling
; global centerDatasetAssembly

; centerDatasetAssembly:



; 	push		ebp
; 	mov			ebp, esp
; 	push		ebx
; 	push		esi
; 	push		edi


; 	mov eax,[ebp+input]    ;dataset
; 	mov ebx,[ebp+input+4]  ;n*4
; 	mov ecx,[ebp+input+4]  ;n*4
; 	mov edx,[ebp+input+12] ;medieDs


; 	xorps xmm4,xmm4
; 	xorps xmm5,xmm5
; 	xorps xmm6,xmm6
; 	xorps xmm7,xmm7


; 	xorps xmm7,xmm7   ;accumulatore finale

; 	mov esi,0


; cicloMedieUnrolling64:

; 			cmp ecx, 64
; 			jl cicloMedie16

; 			movaps xmm0,[eax+esi]	 ;dataset
; 			movaps xmm1,[eax+esi+16]
; 			movaps xmm2,[eax+esi+32]
; 			movaps xmm3,[eax+esi+48]

; 			;evito più che posso i conflitti
; 			addps xmm4,xmm0
; 			addps xmm5,xmm1
; 			addps xmm6,xmm2
; 			addps xmm7,xmm3


; 			add esi,64
; 			sub ecx,64
; 			jmp cicloMedieUnrolling64


; cicloMedie16:

; 			cmp ecx, 16
; 			jl cicloMedieRestoInit
;       		movaps xmm0,[eax+esi]
; 			addps xmm4,xmm0
; 			add esi, 16
; 			sub ecx, 16
; 			jmp cicloMedie16


; cicloMedieRestoInit:

;       xorps xmm1,xmm1 ;riutilizzo xmm1

; cicloMedieResto:

; 			cmp ecx, 4
; 			jl centraMedieFine
; 			movss xmm0,[eax+esi]
; 			movss xmm1,[eax+esi]
; 			add esi, 4
; 			sub ecx, 4
; 			jmp cicloMedieResto


; centraMedieFine:

; 		;a questo punto ho i risultati parziali in xmm4,xmm5,xmm6,xmm7 e xmm1
; 		;li aggiungo in xmm7

; 		addps xmm7,xmm4
; 		addps xmm7,xmm5
; 		addps xmm7,xmm6
; 		addps xmm7,xmm1

; 		;riduzione parallela su xmm7
; 		haddps xmm7,xmm7
; 		haddps xmm7,xmm7

;       	;divido
; 		divss xmm7,[ebp+input+16]

; 		movss [edx],xmm7


; 		;in xmm7 parte bassa abbiamo la media
; 		mov esi,0
; 		shufps xmm7,xmm7, 00000000


; cicloAggDS64:
; 				cmp ebx, 64
; 				jl cicloAggDS16

;  				movaps xmm0,[eax+esi]
;  				movaps xmm1,[eax+esi+16]
;  				movaps xmm2,[eax+esi+32]
;  				movaps xmm3,[eax+esi+48]

; 				subps xmm0,xmm7
; 				subps xmm1,xmm7
; 				subps xmm2,xmm7
; 				subps xmm3,xmm7

; 				movaps [eax+esi],xmm0 ;aggiornamento ds

; 				movaps [eax+esi+16],xmm1 ;aggiornamento ds

; 				movaps [eax+esi+32],xmm2 ;aggiornamento ds

; 				movaps [eax+esi+48],xmm3 ;aggiornamento ds

; 				add esi, 64
; 				sub ebx, 64

; 				jmp cicloAggDS64


;  cicloAggDS16:  ;si sottrae a ogni elemento (scandito per colonna) la media

; 				cmp ebx, 16
; 				jl cicloScalareAggDS

;  				movaps xmm0,[eax+esi]
; 				subps xmm0,xmm7

; 				movaps [eax+esi],xmm0 ;aggiornamento ds
; 				add esi,16
; 				sub ebx, 16

; 				jmp cicloAggDS16

; cicloScalareAggDS:

; 			cmp ebx, 4
; 			jl fineDS

; 			movss xmm0, [eax+esi]
; 			subss xmm0, xmm7

; 			movss [eax+esi], xmm0
; 			add esi, 4
; 			sub ebx, 4

; 			jmp cicloScalareAggDS

; fineDS:

; 	pop	edi
; 	pop	esi
; 	pop	ebx
; 	mov	esp, ebp
; 	pop	ebp
; 	ret






global normalizeVAssembly


normalizeVAssembly:

	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi


	mov eax, [ebp+input] ;v
	mov ebx, [ebp+input+4] ; k*4
	movss xmm7, [ebp+input+8]  ;norma

	sqrtss  xmm7, xmm7
	shufps xmm7,xmm7, 00000000
	mov esi, 0

	cicloUnrollNormalizza64:

						cmp ebx, 64
						jl cicloUnrollNormalizza16

						movaps xmm0, [eax+esi]
						movaps xmm1, [eax+esi+16]
						movaps xmm2, [eax+esi+32]
						movaps xmm3, [eax+esi+48]

						divps xmm0, xmm7
						divps xmm1, xmm7
						divps xmm2, xmm7
						divps xmm3, xmm7

						movaps [eax+esi], xmm0
						movaps [eax+esi+16], xmm1
						movaps [eax+esi+32], xmm2
						movaps [eax+esi+48], xmm3

						add esi, 64
						sub ebx, 64

						jmp cicloUnrollNormalizza64

	cicloUnrollNormalizza16:

						cmp ebx, 16
						jl cicloScalareNormalizza

						movaps xmm0, [eax+esi]
						divps xmm0, xmm7

						movaps [eax+esi], xmm0

						add esi, 16
						sub ebx, 16

						jmp cicloUnrollNormalizza16

	cicloScalareNormalizza:

						; cmp ebx, 4
						; jl fineNormalizza

						; movss xmm0, [eax+esi]
						; divss xmm0, xmm7

						; movss [eax+esi], xmm0

						; add esi, 4
						; sub ebx, 4

						; jmp cicloScalareNormalizza

	fineNormalizza :




	pop	edi
	pop	esi
	pop	ebx
	mov	esp, ebp
	pop	ebp
	ret




global prodottoDsUV2

prodottoDsUV2:



			push		ebp
			mov			ebp, esp
			push		ebx
			push		esi
			push		edi

			mov eax, [ebp+input]    ;dataset
			mov ebx, [ebp+input+4]  ;u
			mov edx, [ebp+input+8]  ;n*4


			mov esi,0
			mov ecx, 0      ;indice su u
			mov edi,edx     ;edi=n
			xorps xmm1,xmm1 ;accumlatore v[0]
			xorps xmm2,xmm2 ;accumlatore v[1]
			xorps xmm3,xmm3 ;accumlatore v[2]
			xorps xmm4,xmm4 ;accumlatore v[3]


cicloPrimaColonnaProdottoDsUV2:

			cmp esi,edx
			jg cicloSecondaColonnaProdottoDsUV2Init

			movaps xmm0, [eax+esi]  ;ds[i..4]
			mulps  xmm0, [ebx+ecx]  ;ds[i..4]*u[i..4]
			addps  xmm1, xmm0		;v[0]


			add esi,16
			add ecx,16
			jmp cicloPrimaColonnaProdottoDsUV2




cicloSecondaColonnaProdottoDsUV2Init:

			sub esi,16
			add edx,edi
			mov ecx,0

cicloSecondaColonnaProdottoDsUV2:

			cmp esi,edx
			jg cicloTerzaColonnaProdottoDsUV2Init

			movaps xmm0, [eax+esi]  ;ds[i..4]
			mulps  xmm0, [ebx+ecx]  ;ds[i..4]*u[i..4]
			addps  xmm2, xmm0		;v[1]


			add esi,16
			add ecx,16
			jmp cicloSecondaColonnaProdottoDsUV2


cicloTerzaColonnaProdottoDsUV2Init:
			sub esi,16
			add edx,edi
			mov ecx,0

cicloTerzaColonnaProdottoDsUV2:

			cmp esi,edx
			jg cicloQuartaColonnaProdottoDsUV2Init

			movaps xmm0, [eax+esi]  ;ds[i..4]
			mulps  xmm0, [ebx+ecx]  ;ds[i..4]*u[i..4]
			addps  xmm3, xmm0		;v[2]


			add esi,16
			add ecx,16
			jmp cicloTerzaColonnaProdottoDsUV2



cicloQuartaColonnaProdottoDsUV2Init:

			sub esi,16
			add edx,edi
			mov ecx,0

cicloQuartaColonnaProdottoDsUV2:

			cmp esi,edx
			jg fineProdottoDsUV2

			movaps xmm0, [eax+esi]  ;ds[i..4]
			mulps  xmm0, [ebx+ecx]  ;ds[i..4]*u[i..4]
			addps  xmm4, xmm0		;v[3]


			add esi,16
			add ecx,16
			jmp cicloQuartaColonnaProdottoDsUV2


fineProdottoDsUV2:




			haddps xmm1,xmm1
			haddps xmm2,xmm2
			haddps xmm3,xmm3
			haddps xmm4,xmm4

			haddps xmm1,xmm1
			haddps xmm2,xmm2
			haddps xmm3,xmm3
			haddps xmm4,xmm4

			xorps xmm7,xmm7



			insertps xmm7, xmm1, 00000000b
			insertps xmm7, xmm2, 00010000b
			insertps xmm7, xmm3, 00100000b
			insertps xmm7, xmm4, 00110000b

			movss xmm5,[ebp+input+16]
			shufps xmm5,xmm5, 00000000
			divps xmm7,xmm5

			mov ecx, [ebp+input+12] ;v
			movaps [ecx], xmm7




			pop	edi
			pop	esi
			pop	ebx
			mov	esp, ebp
			pop	ebp
			ret





global updateU

updateU:


			push		ebp
			mov			ebp, esp
			push		ebx
			push		esi
			push		edi

			mov eax, [ebp+input] ;dataset
			mov ebx, [ebp+input+4] 	;u
			movss xmm6, [ebp+input+8]  ;v[p]
			mov edx, [ebp+input+12] ;n*4
			mov edi, edx ;n*4


			mov esi,0 ;indice i

			shufps xmm6,xmm6, 00000000 ;per moltiplicare v per tutta la colonna replico v


cicloColonneDataset64:


			cmp edi,64
			jl cicloColonneDataset16

			movaps xmm0, [eax+esi]  ;leggi 4 elem da dataset
			movaps xmm1, [eax+esi+16]
			movaps xmm2, [eax+esi+32]
			movaps xmm3, [eax+esi+48]

			mulps xmm0, xmm6   ;effettua prodotto ds*v
			mulps xmm1, xmm6
			mulps xmm2, xmm6
			mulps xmm3, xmm6

			addps xmm0, [ebx+esi]   ;aggiungi a u[i]
			addps xmm1, [ebx+esi+16]
			addps xmm2, [ebx+esi+32]
			addps xmm3, [ebx+esi+48]

			movaps [ebx+esi],xmm0   ;scrivi u[i] in memoria
			movaps [ebx+esi+16],xmm1   ;scrivi u[i] in memoria
			movaps [ebx+esi+32],xmm2   ;scrivi u[i] in memoria
			movaps [ebx+esi+48],xmm3   ;scrivi u[i] in memoria

			add esi,64
			sub edi,64
			jmp cicloColonneDataset64




cicloColonneDataset16:

			cmp edi,16
			jl cicloColonneDatasetResto
			movaps xmm0, [eax+esi]  ;leggi 4 elem da dataset
			mulps  xmm0,xmm6        ;v[p]*dataset[i..4]
			addps xmm0, [ebx+esi]	;leggi 4 elem da u e aggiungi
			movaps [ebx+esi],xmm0   ;scrivi u[i] in memoria

			add esi,16
			sub edi,16
			jmp cicloColonneDataset16


cicloColonneDatasetResto:

			; cmp edi,4
			; jl fineupdateU
			; movss xmm0, [eax+esi]
			; mulss xmm0,xmm6
			; addss xmm0, [ebx+esi]


			; movss [ebx+esi],xmm0

			; add esi,4
			; sub edi,4
			; jmp cicloColonneDatasetResto



fineupdateU:

			pop	edi
			pop	esi
			pop	ebx
			mov	esp, ebp
			pop	ebp
			ret



global resetVectorAssembly

resetVectorAssembly:

			push		ebp
			mov			ebp, esp
			push		ebx
			push		esi
			push		edi


			mov eax, [ebp+input]   ;u
			mov ebx, [ebp+input+4] ;n*4
			mov edx, ebx 	       ;n*4

			xorps xmm0,xmm0
			mov esi,0

cicloAzzeraVettore64:

			cmp edx, 64
			jl cicloAzzeraVettore

			movaps [eax+esi], xmm0
			movaps [eax+esi+16], xmm0
			movaps [eax+esi+32], xmm0
			movaps [eax+esi+48], xmm0

			add esi, 64
			sub edx, 64
			jmp cicloAzzeraVettore64

cicloAzzeraVettore:

			cmp edx,16
			jl cicloRestoAzzeraVettore
			movaps [eax+esi],xmm0
			add esi,16
			sub edx,16
			jmp cicloAzzeraVettore


cicloRestoAzzeraVettore:

			; cmp edx,4
			; jl fineAzzeraVettore
			; movss [eax+esi],xmm0
			; add esi,4
			; sub edx,4
			; jmp cicloRestoAzzeraVettore

fineAzzeraVettore:

			pop	edi
			pop	esi
			pop	ebx
			mov	esp, ebp
			pop	ebp
			ret


global prodottoDsU

prodottoDsU:


			push		ebp
			mov			ebp, esp
			push		ebx
			push		esi
			push		edi


			mov eax, [ebp+input]   ;dataset
			mov ebx, [ebp+input+4] ;u
			mov edx, [ebp+input+8] ;n*4

			mov esi,0
			mov edi,edx

			xorps xmm7,xmm7 ;accumulatore somma
			xorps xmm6,xmm6

cicloProdottoDsU16:

			cmp edi,16
			jl cicloProdottoDsUScalare
			movaps xmm0,[eax+esi]
			mulps  xmm0,[ebx+esi]
			addps  xmm7,xmm0
			add esi,16
			sub edi,16
			jmp cicloProdottoDsU16



cicloProdottoDsUScalare:

			cmp edi,4
			jl fineProdottoDsU
			movaps xmm0,[eax+esi]
			mulss  xmm0,[ebx+esi]
			addss  xmm6,xmm0
			add esi,4
			sub edi,4
			jmp cicloProdottoDsU16


fineProdottoDsU:
			addps xmm7,xmm6
			haddps xmm7,xmm7
			haddps xmm7,xmm7



			getmem 4,1

			movss [eax],xmm7



			pop	edi
			pop	esi
			pop	ebx
			mov	esp, ebp
			pop	ebp
			ret



global updateMatrix

updateMatrix:


            push        ebp
            mov            ebp, esp
            push        ebx
            push        esi
            push        edi

            mov eax, [ebp+input]         ; U+j*input->n
            mov ebx, [ebp+input+4]         ;u
            mov ecx, [ebp+input+8]      ;input->n*4

            mov esi, 0

cicloupdateMatrix64:
			cmp ecx, 64
			jl cicloupdateMatrix16

			movaps xmm0, [ebx+esi]
			movaps xmm1, [ebx+esi+16]
			movaps xmm2, [ebx+esi+32]
			movaps xmm3, [ebx+esi+48]

			movaps [eax+esi], xmm0
			movaps [eax+esi+16], xmm1
			movaps [eax+esi+32], xmm2
			movaps [eax+esi+48], xmm3

			add esi, 64
			sub ecx, 64
			jmp cicloupdateMatrix64


cicloupdateMatrix16:

            cmp ecx, 16
            jl cicloupdateMatrixScalare

            movaps xmm0, [ebx+esi]        ;u[1..4]
            movaps [eax+esi], xmm0

            add esi, 16
            sub ecx, 16
            jmp cicloupdateMatrix16

cicloupdateMatrixScalare:

;             cmp ecx, 4
;             jl fineupdateMatrix

;             movss xmm0, [ebx+esi]
;             movss [eax+esi], xmm0

;             add esi, 4
;             sub ecx, 4
;             jmp cicloupdateMatrixScalare

fineupdateMatrix:

            pop    edi
            pop    esi
            pop    ebx
            mov    esp, ebp
            pop    ebp
            ret


global updateDataset

updateDataset:


			push        ebp
            mov         ebp, esp
            push        ebx
            push        esi
            push        edi

			mov eax, [ebp+input]         ;dataset
            mov ebx, [ebp+input+4]       ;u
			mov ecx, [ebp+input+8]       ;v
            mov edi, [ebp+input+16]      ;k
			mov esi, [ebp+input+20] 	 ;n*k


cicloKupdateDatasetVect:
            mov edx, [ebp+input+12]      ;n*4
			cmp edi,4
			jl fineupdateDataset
			sub edi,4
			movss xmm6,[ecx+edi]    	;leggi ultimo elem di v
			shufps xmm6,xmm6,00000000

cicloNupdateDatasetUnrolling64:

			cmp edx,64
			jl cicloNupdateDatasetVect

			sub edx,64     ;n*4 - 4 
			sub esi,64     ;n*k - 4 

			movaps xmm0, [ebx+edx] 	 ;u
			movaps xmm1, [ebx+edx+16]
			movaps xmm2, [ebx+edx+32]
			movaps xmm3, [ebx+edx+48]


			mulps  xmm0,xmm6	;v*u
			mulps  xmm1,xmm6
			mulps  xmm2,xmm6
			mulps  xmm3,xmm6


			movaps xmm4,[eax+esi]  ;D
			movaps xmm5,[eax+esi+16]
			movaps xmm6,[eax+esi+32]
			movaps xmm7,[eax+esi+48]


			subps xmm4,xmm0
			subps xmm5,xmm1
			subps xmm6,xmm2
			subps xmm7, xmm3

			movaps [eax+esi], xmm4 ;scrivi in D
			movaps [eax+esi+16], xmm5 ;scrivi in D
			movaps [eax+esi+32], xmm6 ;scrivi in D
			movaps [eax+esi+48], xmm7

			movss xmm6,[ecx+edi]    	;leggi ultimo elem di v
			shufps xmm6,xmm6,00000000

			jmp cicloNupdateDatasetUnrolling64

cicloNupdateDatasetVect:

			cmp edx,16
			jl cicloNupdateDatasetResto

			sub edx,16
			sub esi,16



			movaps xmm0,[eax+esi]  	;leggi ultimi 4 elem ds
			movaps xmm1, [ebx+edx] 		;leggi ultimi 4 elem u
			mulps  xmm1,xmm6	   		;u[1..4]*v[p]
			subps xmm0,xmm1 	   		;D-u[1..4]*v[p]

			movaps [eax+esi], xmm0 ;scrivi in D
			jmp cicloNupdateDatasetVect


cicloNupdateDatasetResto:

			cmp edx,4
			jl cicloKupdateDatasetVect

			sub edx,4
			sub esi,4

			movss xmm0,[eax+esi]
			movss xmm1, [ebx+edx]
			mulss  xmm1,xmm6
			subss xmm0,xmm1

			movss [eax+esi], xmm0
			jmp cicloNupdateDatasetResto




fineupdateDataset:


			pop    edi
            pop    esi
            pop    ebx
            mov    esp, ebp
            pop    ebp
            ret

 global queryPointASM

queryPointASM:


			push        ebp
            mov         ebp, esp
            push        ebx
            push        esi
            push        edi

			mov eax, [ebp+input]         ;queryPoint
            mov ebx, [ebp+input+4]       ;QS
			mov ecx, [ebp+input+8]       ;medieDS
            mov edi, [ebp+input+12]      ;k*4


            mov esi, 0

cicloQueryPoint64:

			cmp edi, 64
			jl cicloQueryPoint16

			movaps xmm0, [ebx+esi]
			movaps xmm1, [ebx+esi+16]
			movaps xmm2, [ebx+esi+32]
			movaps xmm3, [ebx+esi+48]

			subps xmm0, [ecx+esi]
			subps xmm1, [ecx+esi+16]
			subps xmm2, [ecx+esi+32]
			subps xmm3, [ecx+esi+48]

			movaps [eax+esi], xmm0
			movaps [eax+esi+16], xmm1
			movaps [eax+esi+32], xmm2
			movaps [eax+esi+48], xmm3

			add esi, 64
			sub edi, 64
			jmp cicloQueryPoint64


cicloQueryPoint16:

            cmp edi, 16
            jl cicloQueryPointResto

            movaps xmm0, [ebx+esi]		;QS
            subps xmm0, [ecx+esi]		; QS-medieDS

            movaps [eax+esi], xmm0

            add esi, 16
            sub edi, 16
            jmp cicloQueryPoint16

cicloQueryPointResto:
			; cmp edi, 4
			; jl fineQP

			; movss xmm0, [ebx+esi]
			; subss xmm0, [ecx+esi]
			; movss [eax+esi], xmm0

			; add esi, 4
			; sub edi, 4
			; jmp cicloQueryPointResto

fineQP:



            pop    edi
            pop    esi
            pop    ebx
            mov    esp, ebp
            pop    ebp
            ret




global queryPointASMNoPCA

queryPointASMNoPCA:


			push        ebp
            mov         ebp, esp
            push        ebx
            push        esi
            push        edi

						mov eax, [ebp+input]         ;queryPoint
            mov ebx, [ebp+input+4]       ;QS
            mov edi, [ebp+input+8]      ;k*4


            mov esi, 0

cicloQueryPoint64NoPCA:

			cmp edi, 64
			jl cicloQueryPoint16NoPCA

			movaps xmm0, [ebx+esi]
			movaps xmm1, [ebx+esi+16]
			movaps xmm2, [ebx+esi+32]
			movaps xmm3, [ebx+esi+48]

			movaps [eax+esi], xmm0
			movaps [eax+esi+16], xmm1
			movaps [eax+esi+32], xmm2
			movaps [eax+esi+48], xmm3

			add esi, 64
			sub edi, 64

			jmp cicloQueryPoint64NoPCA


cicloQueryPoint16NoPCA:

            cmp edi, 16
            jl cicloQueryPointRestoNoPCA

            movaps xmm0, [ebx+esi]		;QS

            movaps [eax+esi], xmm0

            add esi, 16
            sub edi, 16
            jmp cicloQueryPoint16NoPCA

cicloQueryPointRestoNoPCA:

			; cmp edi, 4
			; jl fineQPNoPCA

			; movss xmm0, [ebx+esi]
			; movss [eax+esi], xmm0

			; add esi, 4
			; sub edi, 4
			; jmp cicloQueryPointRestoNoPCA

fineQPNoPCA:


            pop    edi
            pop    esi
            pop    ebx
            mov    esp, ebp
            pop    ebp
            ret

global productQPV

productQPV:


			push        ebp
            mov         ebp, esp
            push        ebx
            push        esi
            push        edi

			mov eax, [ebp+input]         ;queryPoint
            mov ebx, [ebp+input+4]       ;V+i*input->k
            mov edi, [ebp+input+8]      ;k*4
			mov ecx, [ebp+input+12]

            mov esi,0
            xorps xmm7, xmm7
            xorps xmm6, xmm6


cicloproductQPV64:

			cmp edi, 64
			jl cicloproductQPV16

			movaps xmm0, [eax+esi]
			movaps xmm1, [eax+esi+16]
			movaps xmm2, [eax+esi+32]
			movaps xmm3, [eax+esi+48]

			mulps xmm0, [ebx+esi]
			mulps xmm1, [ebx+esi+16]
			mulps xmm2, [ebx+esi+32]
			mulps xmm3, [ebx+esi+48]

			addps xmm7, xmm0
			addps xmm7, xmm1
			addps xmm7, xmm2
			addps xmm7, xmm3

			add esi, 64
			sub edi, 64
			jmp cicloproductQPV64

cicloproductQPV16:

           cmp edi, 16
           jl cicloproductQPVScalare

           movaps xmm0, [eax+esi]
           mulps xmm0, [ebx+esi]
           addps xmm7, xmm0

           sub edi, 16
           add esi, 16

           jmp cicloproductQPV16


cicloproductQPVScalare:

			cmp edi, 4
			jl fineproductQPV

			movss xmm1, [eax+esi]
			mulss xmm1, [ebx+esi]

			addps xmm6, xmm1

			sub edi, 4
			add esi, 4

			jmp cicloproductQPVScalare

fineproductQPV:

			addps xmm7, xmm6
			haddps xmm7, xmm7
			haddps xmm7, xmm7


			movss [ecx], xmm7


            pop    edi
            pop    esi
            pop    ebx
            mov    esp, ebp
            pop    ebp
            ret
