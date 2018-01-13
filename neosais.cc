#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>

#ifndef IMP_SIGMA
#define IMP_SIGMA 256
#endif

#define MAXN UINT32_MAX - 1

#define T_char unsigned char
#define T_uint uint32_t
#define T_iter int64_t
#define T_size int64_t
#define T_type bool
#define T_flag bool

		inline bool is_lms(const T_type *T, const T_iter x){return x && ~x && T[x] && !T[x-1];}

		inline bool is_legal(const T_uint x){return x && ~x;}


	template <typename T_data>
	inline bool equal_str(
		const T_data *S, const T_type *T,
		T_uint x, T_uint y, const T_size n)
	{
		do{if(S[x++] != S[y++]) return 0;}
		while(!is_lms(T, x) && !is_lms(T, y) && x < n && y < n);
		if(x == n || y == n) return 0;
		return S[x] == S[y];
	}

template <typename T_data>
int neo_sais_core(
	T_data *S, T_uint *A,
	const T_size n, const T_size k,
	T_type *T, T_uint *P)
{
	T_uint *B, *C, *RS, *RA;
	T_iter i, j, lx = -1;
	T_size m = 0, nm;

	T[n-1] = 0;
	for(i = n-2; i >= 0; i--) T[i] = S[i] < S[i+1] || (S[i] == S[i+1] && T[i+1]);
	for(i = 1; i < n; i++) if(T[i] && !T[i-1]) P[m++] = i;

	B = P + m;
	memset(B, 0, sizeof(T_uint) * k);
	for(i = 0; i < n; i++) B[S[i]]++;
	for(i = 1; i < k; i++) B[i] += B[i-1];

	C = B + k;
	memcpy(C, B, sizeof(T_uint) * k);

	memset(A, T_uint(-1), sizeof(T_uint) * n);
	for(i = 0; i < m; i++) A[--C[S[P[i]]]] = P[i];

	C[0] = 0; memcpy(C + 1, B, sizeof(T_uint) * (k - 1));
	A[C[S[n-1]]++] = n-1;
	for(i = 0; i < n; i++) if(is_legal(A[i]) && !T[A[i]-1]) A[C[S[A[i]-1]]++] = A[i]-1;
	memcpy(C, B, sizeof(T_uint) * k);
	for(i = n-1; i >= 0; i--) if(is_legal(A[i]) && T[A[i]-1]) A[--C[S[A[i]-1]]] = A[i]-1;

	if(m){
	RS = C + k; RA = RS + m; j = nm = 0;
	memset(RS, T_uint(-1), sizeof(T_uint) * n >> 1);
	while(!is_lms(T, A[j])) j++;
	lx = A[j]; RS[A[j] >> 1] = nm;
	for(i = j + 1; i < n; i++)
	if(is_lms(T, A[i])){
		if(!equal_str(S, T, A[i], lx, n)) nm++;
		lx = A[i]; RS[A[i] >> 1] = nm;
	}nm++;
	for(i = j = 0; i <= n >> 1; i++) if(~RS[i]) RS[j++] = RS[i];

 	if(nm < m) neo_sais_core(RS, RA, m, nm, T+n, RA+m);
	else for(i = 0; i < m; i++) RA[RS[i]] = i;

	memset(A, T_uint(-1), sizeof(T_uint) * n);
	memcpy(C, B, sizeof(T_uint) * k);
	for(i = m-1; i >= 0; i--) A[--C[S[P[RA[i]]]]] = P[RA[i]];

	C[0] = 0; memcpy(C + 1, B, sizeof(T_uint) * (k - 1));
	A[C[S[n-1]]++] = n-1;
	for(i = 0; i < n; i++) if(is_legal(A[i]) && !T[A[i]-1]) A[C[S[A[i]-1]]++] = A[i]-1;

	memcpy(C, B, sizeof(T_uint) * k);
	for(i = n-1; i >= 0; i--) if(is_legal(A[i]) && T[A[i]-1]) A[--C[S[A[i]-1]]] = A[i]-1;
	}
	return 0;
}

template <typename T_data>
int neo_sais_uchar(
	T_data *S, T_uint *A,
	const T_size n, const T_size k = IMP_SIGMA)
{
	if(S == NULL) return -1;
	if(A == NULL) return -2;
	if(n < 1 || n > MAXN) return -3;
	if(k < 1 || k > IMP_SIGMA) return -4;
	T_uint *P; T_type *T;
	if(NULL == (P =(T_uint *)malloc(sizeof(T_uint) * std::max(n, k) << 2))) return -5;
	if(NULL == (T =(T_type *)malloc(sizeof(T_type) * n << 1))) return -6;

	if(n == 1) {A[0] = 0; return 0;}

	neo_sais_core(S, A, n, k, T, P);
	free(P); free(T);
	return 0;
}
