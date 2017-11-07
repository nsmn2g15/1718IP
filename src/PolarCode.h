/*    
PolarCode provides functions for polar code training, encoding and decoder.
Copyright (C) 2010  Robert G. Maunder

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

The GNU General Public License can be seen at http://www.gnu.org/licenses/.
*/


/*! \mainpage PolarCode 1.1.0
 * PolarCode provides functions for polar code training, encoding and decoding. Like turbo and LDPC codes, polar codes facilitate near-capacity operation. However, unlike turbo and LDPC codes, polar codes do not require an iterative decoder. You can find out more about polar codes in <A href="http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5075875">Erdal Arikan’s paper</A>. Tested with <A href="http://itpp.sourceforge.net/current/index.html">IT++</A> 4.0.2 and g++ 4.2.1.
 * 
 * Get started in the <A href="group__polar.html">PolarCode page</A>.
 * 
 * Copyright &copy; 2010 Robert G. Maunder. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the <A href="http://www.gnu.org/licenses/">GNU General Public License</A> for more details.
 */

/// \file PolarCode.h
/// \author Robert G Maunder
/// \date 03/08/2010
/// \version 1.1.0
/// \brief PolarCode provides functions for polar code training, encoding and decoder.


#ifndef _POLARCODE_H
#define _POLARCODE_H

#include <iostream>
#include "itpp/itbase.h"
#include "robprob.h"

using namespace std;
using namespace itpp;
using namespace RobProb;

/// \defgroup polar PolarCode


/// \ingroup polar
/// Class for polar encoding and decoding.
class PolarCode
{
	public:	
		/// This constructor allows the parameters of the polar code to be set.
		/// \param new_N The encoded block length \f$N\f$.
		/// \param new_K The uncoded block length \f$K\f$.
		/// \param new_A The information set \f$\mathcal{A}\f$.
		/// \param new_u_A_c The frozen bits \f$u_{\mathcal{A}^c}\f$. If this is not provided, then the frozen bits will be randomised.
		inline PolarCode(int new_N, int new_K = 0, const ivec& new_A = "", const bvec& new_u_A_c = "")
		{
			set(new_N, new_K, new_A, new_u_A_c);
		}
				
		/// This constructor allows the parameters of the polar code to be read from a parser. For example, "N=8 K=5 A=[3,5,6,7,8] u_A_c=[0,0,0]". If the frozen bits cannot be read from the parser, they will be randomised.
		/// \param parser The parser.
		inline PolarCode(Parser& parser)
		{
			set(parser);
		}
		
		
		/// This function allows the parameters of the polar code to be set.
		/// \param new_N The encoded block length \f$N\f$.
		/// \param new_K The information block length \f$K\f$.
		/// \param new_A The information set \f$\mathcal{A}\f$.
		/// \param new_u_A_c The frozen bits \f$u_{\mathcal{A}^c}\f$. If this is not provided, then the frozen bits will be randomised.
		void set(int new_N, int new_K = 0, const ivec& new_A = "", const bvec& new_u_A_c = "");	
		
		/// This function allows the parameters of the polar code to be read from a parser. For example, "N=8 K=5 A=[3,5,6,7,8] u_A_c=[0,0,0]". If the frozen bits cannot be read from the parser, they will be randomised.
		/// \param parser The parser.
		void set(Parser& parser);		
		
		/// This function sets the frozen bits \f$u_{\mathcal{A}^c}\f$.
		/// \param new_u_A_c The frozen bits \f$u_{\mathcal{A}^c}\f$.
		void set_u_A_c(const bvec& new_u_A_c);
		
		/// This function randomises the frozen bits \f$u_{\mathcal{A}^c}\f$.
		inline void randomise_u_A_c(void)
		{
			u_A_c = randb(N-K);
		}
		
		/// This function provides the encoded block length \f$N\f$.
		inline int get_N(void) const
		{
			return N;
		}
		
		/// This function provides the information block length \f$K\f$.
		inline int get_K(void) const
		{
			return K;
		}
		
		/// This function provides the information set \f$\mathcal{A}\f$.
		inline const ivec& get_A(void) const
		{
			return A;
		}
		
		/// This function provides the frozen set \f$\mathcal{A}^c\f$.
		inline const ivec& get_A_c(void) const
		{
			return A_c;
		}
		
		/// This function provides the frozen bits \f$u_{\mathcal{A}^c}\f$.
		inline const bvec& get_u_A_c(void) const
		{
			return u_A_c;
		}
				
		/// This function provides the code rate \f$K/N\f$.
		inline double get_code_rate(void) const
		{
			return static_cast<double>(K)/static_cast<double>(N);
		}
		
		/// This function provides the generator matrix \f$G_N\f$.
		/// \return The generator matrix \f$G_N\f$.
		inline const bmat& get_G_N(void) const
		{
			return G_N;
		}
		
		/// This function performs the polar encoding operation \f$x_1^N = u_\mathcal{A}G_N(\mathcal{A}) \oplus u_{\mathcal{A}^c}G_N(\mathcal{A}^c)\f$.
		/// \param u_A The information bits \f$u_\mathcal{A}\f$.
		/// \return The encoded bits \f$x_1^N\f$.
		inline bvec encode(const bvec& u_A = "") const
		{
			bvec x;
			encode(u_A, x);
			return x;
		}
		
		/// This function performs the polar encoding operation \f$x_1^N = u_\mathcal{A}G_N(\mathcal{A}) \oplus u_{\mathcal{A}^c}G_N(\mathcal{A}^c)\f$.
		/// \param u_A The information bits \f$u_\mathcal{A}\f$.
		/// \param x The encoded bits \f$x_1^N\f$.
		void encode(const bvec& u_A, bvec& x) const;
		
		void decode_slow(const plr_frame& L_1, plr_frame& L_N, bvec& u_hat) const;
	
		inline void decode_slow(const llr_frame& LL_1, llr_frame& LL_N, bvec& u_hat) const
		{
			plr_frame L_N;			
			decode_slow(exp(LL_1), L_N, u_hat);
			LL_N = log(L_N);
		}
		
		/// This function performs polar decoding.
		/// \param LL_1 The logarithmic likelihood ratios \f$ \{\ln L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		/// \return The logarithmic likelihood ratios \f$ \{\ln L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):i \in \mathcal{A}\} \f$
		inline llr_frame decode(const llr_frame& LL_1, int& acs_operations = dummy_int) const
		{
			llr_frame LL_A_N;
			bvec u_hat;
			decode(LL_1, LL_A_N, u_hat, acs_operations);
			return LL_A_N;
		}
		
		/// This function performs polar decoding.
		/// \param LL_1 The logarithmic likelihood ratios \f$ \{\ln L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param LL_A_N The logarithmic likelihood ratios \f$ \{\ln L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):i \in \mathcal{A}\} \f$
		/// \param u_A_hat The decoded bits \f$ \hat{u}_\mathcal{A}\f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		void decode(const llr_frame& LL_1, llr_frame& LL_A_N, bvec& u_A_hat = dummy_bvec, int& acs_operations = dummy_int) const;
		
		/// This function performs polar decoding.
		/// \param LL_1 The logarithmic likelihood ratios \f$ \{\ln L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param LL_N The logarithmic likelihood ratios \f$ \{\ln L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):1 \leq i \leq N\} \f$
		/// \param u_hat The decoded bits \f$ \hat{u}_1^N\f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		void decode_frozen(const llr_frame& LL_1, llr_frame& LL_N, bvec& u_hat = dummy_bvec, int& acs_operations = dummy_int) const;
		
		/// This function performs polar decoding.
		/// \param LL_1 The logarithmic likelihood ratios \f$ \{\ln L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param u_hat The decoded bits \f$ \hat{u}_1^N\f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		/// \return The logarithmic likelihood ratios \f$ \{\ln L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):1 \leq i \leq N\} \f$
		inline llr_frame decode_frozen(const llr_frame& LL_1, bvec& u_hat = dummy_bvec, int& acs_operations = dummy_int) const
		{
			llr_frame LL_N;
			decode_frozen(LL_1, LL_N, u_hat, acs_operations);
			return LL_N;	
		}
		
		
		
		
		
		
		
		/// This function performs polar decoding.
		/// \param L_1 The likelihood ratios \f$ \{L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		/// \return The likelihood ratios \f$ \{L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):i \in \mathcal{A}\} \f$
		inline plr_frame decode(const plr_frame& L_1, int& acs_operations = dummy_int) const
		{
			plr_frame L_A_N;
			bvec u_hat;
			decode(L_1, L_A_N, u_hat, acs_operations);
			return L_A_N;
		}
				
		/// This function performs polar decoding.
		/// \param L_1 The likelihood ratios \f$ \{L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param L_A_N The likelihood ratios \f$ \{L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):i \in \mathcal{A}\} \f$
		/// \param u_A_hat The decoded bits \f$ \hat{u}_\mathcal{A}\f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		void decode(const plr_frame& L_1, plr_frame& L_A_N, bvec& u_A_hat = dummy_bvec, int& acs_operations = dummy_int) const;
				
		/// This function performs polar decoding.
		/// \param L_1 The likelihood ratios \f$ \{L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param L_N The likelihood ratios \f$ \{L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):1 \leq i \leq N\} \f$
		/// \param u_hat The decoded bits \f$ \hat{u}_1^N\f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		void decode_frozen(const plr_frame& L_1, plr_frame& L_N, bvec& u_hat = dummy_bvec, int& acs_operations = dummy_int) const;
	
		
		/// This function performs polar decoding.
		/// \param L_1 The likelihood ratios \f$ \{L_1^{(1)}(y_i):1 \leq i \leq N\} \f$
		/// \param u_hat The decoded bits \f$ \hat{u}_1^N\f$
		/// \param acs_operations The number of add, compare and select operations performed during decoding.
		/// \return The likelihood ratios \f$ \{L_N^{(i)}(y_1^N,\hat{u}_1^{i-1}):1 \leq i \leq N\} \f$
		inline plr_frame decode_frozen(const plr_frame& L_1, bvec& u_hat = dummy_bvec, int& acs_operations = dummy_int) const
		{
			plr_frame L_N;
			decode_frozen(L_1, L_N, u_hat, acs_operations);
			return L_N;	
		}
		
		
		
		
	private:
		plr L_Nover2_slow(const plr_frame& L_1, const bvec& u_hat) const;		
		
		
		plr L_Nover2(const plr_frame& L_1, const bvec& u_hat, Mat<plr>& Ls, bmat& done, int a, int c=0) const;	
		llr L_Nover2(const llr_frame& L_1, const bvec& u_hat, Mat<llr>& Ls, bmat& done, int a, int c=0) const;	
		
		/// The encoded block length \f$N\f$.
		int N;
		
		/// The information block length \f$K\f$.
		int K;
		
		/// The information set \f$\mathcal{A}\f$.
		ivec A;
		
		/// The frozen set \f$\mathcal{A}^c\f$.
		ivec A_c;
				
		/// The frozen bits \f$u_{\mathcal{A}^c}\f$.
		bvec u_A_c;
		
		/// \f$n = \log_2(N)\f$
		int n;
		
		/// The generator matrix \f$G_N\f$.
		bmat G_N;
		
		
		static int dummy_int;
		
		static bvec dummy_bvec;

	
};

/// \ingroup misc
/// Output a subvector containing the first, third, fifth, etc elements from the input vector
/// \param input The input vector.
/// \return The output vector.
template<typename Num_T> Vec<Num_T> odd(const Vec<Num_T>& input)
{
	Vec<Num_T> output(ceil_i(static_cast<double>(input.length())/2.0));
	
	for(int index = 0; index < output.length(); index++)
	{
		output(index) = input(2*index);
	}
	
	return output;
}

/// \ingroup misc
/// Output a subvector containing the second, fourth, sixth, etc elements from the input vector
/// \param input The input vector.
/// \return The output vector.
template<typename Num_T> Vec<Num_T> even(const Vec<Num_T>& input)
{
	Vec<Num_T> output(floor_i(static_cast<double>(input.length())/2.0));
	
	for(int index = 0; index < output.length(); index++)
	{
		output(index) = input(2*index+1);
	}
	
	return output;
}

#endif
