/*    
PolarCode provides functions for polar code training, encoding and decoder.
Copyright (C) 2010  Robert G. Maunder

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

The GNU General Public License can be seen at http://www.gnu.org/licenses/.
*/

#include "PolarCode.h"

int PolarCode::dummy_int;
bvec PolarCode::dummy_bvec;

void PolarCode::encode(const bvec& u_A, bvec& x) const
{
	string function_name = "void PolarCode::encode(const bvec& u_A, bvec& x) const";
	
	if(u_A.length() != K)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "u_A.length(){" << u_A.length() << "} != K{" << K << "}" << endl;
		exit(1);					
	}
		
	// Perform polar encoding using Equation 9 from the polar code paper.	
	if(A.length() > 0)
	{
		x = (u_A.T()*G_N.get_rows(A-1)+u_A_c.T()*G_N.get_rows(A_c-1)).get_row(0);
	}
	else
	{
		x = (u_A_c.T()*G_N.get_rows(A_c-1)).get_row(0);
	}
	
}

void PolarCode::set(Parser& parser)
{
	string function_name = "void PolarCode::set(Parser& parser)";
		
	int new_N;
	if(parser.exist("N"))
	{
		parser.get(new_N, "N");
	}
	else
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "Cannot find N parameter" << endl;
		exit(1);
	}
		
	int new_K = 0;
	if(parser.exist("K"))
	{
		parser.get(new_K, "K");
	}
	
	ivec new_A = "";
	if(parser.exist("A"))
	{
		parser.get(new_A, "A");
	}
	
	bvec new_u_A_c = "";
	if(parser.exist("u_A_c"))
	{
		parser.get(new_u_A_c, "u_A_c");
	}
		
	set(new_N, new_K, new_A, new_u_A_c);
}



void PolarCode::set(int new_N, int new_K, const ivec& new_A, const bvec& new_u_A_c)
{
	string function_name = "void PolarCode::set(int new_N, int new_K, const ivec& new_A, const bvec& new_u_A_c)";
	
	N = new_N;	
	K = new_K;	
	A = new_A;
	if(A.length() > 1)
	{	
		sort(A);
	}
	if(new_u_A_c.length() == 0)
	{
		randomise_u_A_c();
	}
	else
	{
		u_A_c = new_u_A_c;
	}
		
	if(N < 2)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "N{" << N << "} is not a power of 2" << endl;
		exit(1);			
	}
	
	n = levels2bits(N);
	
	if(pow2i(n) != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "N{" << N << "} is not a power of 2" << endl;
		exit(1);			
	}
	
	if(K < 0)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "K{" << K << "} < 0" << endl;
		exit(1);			
	}
	
	if(K > N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "K{" << K << "} > N{" << N << "}" << endl;
		exit(1);			
	}
	
	if(A.length() != K)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "A.length(){" << A.length() << "} != K{" << K << "}" << endl;
		exit(1);			
	}
	
	if(A.length() > 0)
	{ 
		if(min(A) < 1)
		{
			cerr << "Error!" << endl;
			cerr << function_name << endl;
			cerr << "min(A){" << min(A) << "} < 1" << endl;
			exit(1);			
		}
		
		if(max(A) > N)
		{
			cerr << "Error!" << endl;
			cerr << function_name << endl;
			cerr << "max(A){" << min(A) << "} < N{" << N << "}" << endl;
			exit(1);			
		}
	}
	
	if(A.length() > 1 && min(A(1,A.length()-1) - A(0,A.length()-2)) == 0)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "The elements in A are not unique" << endl;
		exit(1);			
	}
	
	if(u_A_c.length() + K != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "u_A_c.length(){" << u_A_c.length() << "} + K{" << K << "} != N{" << N << "}" << endl;
		exit(1);			
	}
	
	
	A_c.set_size(u_A_c.length());
	
	int A_c_index = 0;
	int A_index = 0;	
	for(int number = 1; number <= N; number++)
	{
		if(A.length() > 0 && A(A_index) == number)
		{
			A_index++;
		}
		else
		{
			A_c(A_c_index) = number;
			A_c_index++;
		}
	}
	
	
	// Set G_N using Equation 72 from the polar code paper.
	G_N.set_size(N, N);
	for(int row_index = 0; row_index < N; row_index++)
	{
		bvec b = dec2bin(n,row_index);
		
		for(int column_index = 0; column_index < N; column_index++)
		{			
			bvec b_dash = dec2bin(n,column_index);
			
			bin value = 1;			
			for(int i = 0; i < n; i++)
			{
				value *= bin(1) + b_dash(i) + b(n-i-1)*b_dash(i);
			}
			G_N(row_index, column_index) = value;
		}
	}
}

void PolarCode::set_u_A_c(const bvec& new_u_A_c)
{
	string function_name = "void PolarCode::set_u_A_c(const bvec& new_u_A_c)";
	
	u_A_c = new_u_A_c;
	
	if(u_A_c.length() + K != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "u_A_c.length(){" << u_A_c.length() << "} + K{" << K << "} != N{" << N << "}" << endl;
		exit(1);			
	}
}

// This was my first attempt at the decoding algorithm. It performs some calculations more than once and so is slower than the optimised function.
/*
void PolarCode::decode_slow(const plr_frame& L_1, plr_frame& L_N, bvec& u_hat) const
{
	string function_name = "void PolarCode::decode_slow(const plr_frame& L_1, plr_frame& L_N, bvec& u_hat) const";
	
	if(L_1.length() != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "L_1.length(){" << L_1.length() << "} != N{" << N << "}" << endl;
		exit(1);			
	}
	
	L_N.set_size(N);
	u_hat.set_size(N);
	
	int A_c_index = 0;
	
	for(int i = 0; i < N; i++)
	{
		if(i > 0)
		{
			L_N(i) = L_Nover2_slow(L_1, u_hat(0,i-1));
		}
		else
		{
			L_N(i) = L_Nover2_slow(L_1, "");
		}
		
		if(A_c_index < u_A_c.size() && i == A_c(A_c_index)-1 )
		{
			u_hat(i) = u_A_c(A_c_index);
			A_c_index++;
		}
		else
		{
			u_hat(i) = hard(L_N(i));
		}
		
	}

}

plr PolarCode::L_Nover2_slow(const plr_frame& L_1, const bvec& u_hat) const
{	
	if(L_1.length() == 1)
	{
		return L_1(0);
	}
	else
	{
		if(u_hat.length() % 2 == 0) // use equation 74
		{
			return ( L_Nover2_slow(L_1(0,L_1.length()/2-1),odd(u_hat)+even(u_hat)) * L_Nover2_slow(L_1(L_1.length()/2,L_1.length()-1),even(u_hat)) + 1 ) / ( L_Nover2_slow(L_1(0,L_1.length()/2-1),odd(u_hat)+even(u_hat)) + L_Nover2_slow(L_1(L_1.length()/2,L_1.length()-1),even(u_hat)) );
			
		}
		else // use equation 75
		{
			if(u_hat(u_hat.length()-1))
			{
				if(u_hat.length() < 2)
				{
					return L_Nover2_slow(L_1(L_1.length()/2,L_1.length()-1),"") / L_Nover2_slow(L_1(0,L_1.length()/2-1),"");
				}
				else
				{
				
					return L_Nover2_slow(L_1(L_1.length()/2,L_1.length()-1),even(u_hat(0,length(u_hat)-2))) / L_Nover2_slow(L_1(0,L_1.length()/2-1),odd(u_hat(0,length(u_hat)-2))+even(u_hat(0,length(u_hat)-2)));
				}
			}
			else
			{
				if(u_hat.length() < 2)
				{
					return L_Nover2_slow(L_1(L_1.length()/2,L_1.length()-1),"") * L_Nover2_slow(L_1(0,L_1.length()/2-1),"");
				}
				else
				{
					return L_Nover2_slow(L_1(L_1.length()/2,L_1.length()-1),even(u_hat(0,length(u_hat)-2))) * L_Nover2_slow(L_1(0,L_1.length()/2-1),odd(u_hat(0,length(u_hat)-2))+even(u_hat(0,length(u_hat)-2)));
				}
			}
		}
	}
}
*/

void PolarCode::decode_frozen(const plr_frame& L_1, plr_frame& L_N, bvec& u_hat, int& acs_operations) const
{	
	string function_name = "void PolarCode::decode(const plr_frame& L_1, plr_frame& L_N, bvec& u_hat, int& acs_operations) const";
	
	if(L_1.length() != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "L_1.length(){" << L_1.length() << "} != N{" << N << "}" << endl;
		exit(1);			
	}
	
	acs_count = 0;
	
	L_N.set_size(N);
	u_hat.set_size(N);
	
	int A_c_index = 0;
	
	Mat<plr> Ls(N,n);
	bmat done(N,n);
	done.zeros();
	
	for(int i = 0; i < N; i++)
	{
		if(i > 0)
		{
			L_N(i) = L_Nover2(L_1, u_hat(0,i-1), Ls, done, i);
		}
		else
		{
			L_N(i) = L_Nover2(L_1, "", Ls, done, i);
		}
		
		if(A_c_index < u_A_c.size() && i == A_c(A_c_index)-1 )
		{
			u_hat(i) = u_A_c(A_c_index);
			A_c_index++;
		}
		else
		{
			u_hat(i) = hard(L_N(i));
		}	
	}

	acs_operations = acs_count;
}

void PolarCode::decode(const plr_frame& L_1, plr_frame& L_A_N, bvec& u_A_hat, int& acs_operations) const
{

	
	string function_name = "void PolarCode::decode(const plr_frame& L_1, plr_frame& L_A_N, bvec& u_A_hat, int& acs_operations) const";
	
	if(L_1.length() != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "L_1.length(){" << L_1.length() << "} != N{" << N << "}" << endl;
		exit(1);			
	}
	
	acs_count = 0;
	
	plr_frame L_N(N);
	bvec u_hat(N);
	
	int A_c_index = 0;
	
	Mat<plr> Ls(N,n);
	bmat done(N,n);
	done.zeros();
	
	for(int i = 0; i < N; i++)
	{
		
		if(A_c_index < u_A_c.size() && i == A_c(A_c_index)-1 )
		{
			u_hat(i) = u_A_c(A_c_index);
			A_c_index++;
		}
		else
		{
			if(i > 0)
			{
				L_N(i) = L_Nover2(L_1, u_hat(0,i-1), Ls, done, i);
			}
			else
			{
				L_N(i) = L_Nover2(L_1, "", Ls, done, i);
			}
			
			u_hat(i) = hard(L_N(i));
		}	
	}
	
	L_A_N = L_N(A-1);
	u_A_hat = u_hat(A-1);
	
	acs_operations = acs_count;
}

plr PolarCode::L_Nover2(const plr_frame& L_1, const bvec& u_hat, Mat<plr>& Ls, bmat& done, int a, int c) const
{	
//	cout << a+1 << '\t' << length(L_1) << '\t' << c+1 << '\t' << N*a/L_1.length() + c/L_1.length() << '\t' << levels2bits(L_1.length())-1 << endl;
	
	if(L_1.length() == 1)
	{
		return L_1(0);
	}
	else
	{
		int x = N*a/L_1.length() + c/L_1.length();
		int y = levels2bits(L_1.length())-1;
		
		if(done(x,y) == 0)
		{
			if(u_hat.length() % 2 == 0) // use equation 74
			{
				plr left = L_Nover2(L_1(0,L_1.length()/2-1),odd(u_hat)+even(u_hat), Ls, done, a/2, c);
				plr right = L_Nover2(L_1(L_1.length()/2,L_1.length()-1),even(u_hat), Ls, done, a/2, c+L_1.length()/2);
				Ls(x,y) = (left*right+1)/(left+right);
				
			}
			else // use equation 75
			{
				bvec left_b_hat = "";
				bvec right_b_hat = "";
	
				if(u_hat.length() >= 2)
				{
					left_b_hat = odd(u_hat(0,length(u_hat)-2))+even(u_hat(0,length(u_hat)-2));
					right_b_hat = even(u_hat(0,length(u_hat)-2));
				}
				
				plr left = L_Nover2(L_1(0,L_1.length()/2-1),left_b_hat, Ls, done, (a-1)/2, c);
				plr right = L_Nover2(L_1(L_1.length()/2,L_1.length()-1),right_b_hat, Ls, done, (a-1)/2, c+L_1.length()/2);
				
				if(u_hat(u_hat.length()-1))
				{
	
					Ls(x,y) = right/left;
				}
				else
				{
					Ls(x,y) = right*left;
				}
			}
			done(x,y) = 1;
		}
			
		return Ls(x,y);
	}
}

void PolarCode::decode_frozen(const llr_frame& L_1, llr_frame& L_N, bvec& u_hat, int& acs_operations) const
{	
	string function_name = "void PolarCode::decode(const llr_frame& L_1, llr_frame& L_N, bvec& u_hat, int& acs_operations) const";
	
	if(L_1.length() != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "L_1.length(){" << L_1.length() << "} != N{" << N << "}" << endl;
		exit(1);			
	}
	
	acs_count = 0;
	
	L_N.set_size(N);
	u_hat.set_size(N);
	
	int A_c_index = 0;
	
	Mat<llr> Ls(N,n);
	bmat done(N,n);
	done.zeros();
	
	for(int i = 0; i < N; i++)
	{
		if(i > 0)
		{
			L_N(i) = L_Nover2(L_1, u_hat(0,i-1), Ls, done, i);
		}
		else
		{
			L_N(i) = L_Nover2(L_1, "", Ls, done, i);
		}
		
		if(A_c_index < u_A_c.size() && i == A_c(A_c_index)-1 )
		{
			u_hat(i) = u_A_c(A_c_index);
			A_c_index++;
		}
		else
		{
			u_hat(i) = hard(L_N(i));
		}	
	}

	acs_operations = acs_count;
}

void PolarCode::decode(const llr_frame& L_1, llr_frame& L_A_N, bvec& u_A_hat, int& acs_operations) const
{

	
	string function_name = "void PolarCode::decode(const llr_frame& L_1, llr_frame& L_A_N, bvec& u_A_hat, int& acs_operations) const";
	
	if(L_1.length() != N)
	{
		cerr << "Error!" << endl;
		cerr << function_name << endl;
		cerr << "L_1.length(){" << L_1.length() << "} != N{" << N << "}" << endl;
		exit(1);			
	}
	
	acs_count = 0;
	
	llr_frame L_N(N);
	bvec u_hat(N);
	
	int A_c_index = 0;
	
	Mat<llr> Ls(N,n);
	bmat done(N,n);
	done.zeros();
	
	for(int i = 0; i < N; i++)
	{
		
		if(A_c_index < u_A_c.size() && i == A_c(A_c_index)-1 )
		{
			u_hat(i) = u_A_c(A_c_index);
			A_c_index++;
		}
		else
		{
			if(i > 0)
			{
				L_N(i) = L_Nover2(L_1, u_hat(0,i-1), Ls, done, i);
			}
			else
			{
				L_N(i) = L_Nover2(L_1, "", Ls, done, i);
			}
			
			u_hat(i) = hard(L_N(i));
		}	
	}
	
	L_A_N = L_N(A-1);
	u_A_hat = u_hat(A-1);
	
	acs_operations = acs_count;
}



llr PolarCode::L_Nover2(const llr_frame& L_1, const bvec& u_hat, Mat<llr>& Ls, bmat& done, int a, int c) const
{	
//	cout << a+1 << '\t' << length(L_1) << '\t' << c+1 << '\t' << N*a/L_1.length() + c/L_1.length() << '\t' << levels2bits(L_1.length())-1 << endl;
	
	if(L_1.length() == 1)
	{
		return L_1(0);
	}
	else
	{
		int x = N*a/L_1.length() + c/L_1.length();
		int y = levels2bits(L_1.length())-1;
		
		if(done(x,y) == 0)
		{
			if(u_hat.length() % 2 == 0) // use equation 74
			{
				llr left = L_Nover2(L_1(0,L_1.length()/2-1),odd(u_hat)+even(u_hat), Ls, done, a/2, c);
				llr right = L_Nover2(L_1(L_1.length()/2,L_1.length()-1),even(u_hat), Ls, done, a/2, c+L_1.length()/2);
				Ls(x,y) = left^right;
				
			}
			else // use equation 75
			{
				bvec left_b_hat = "";
				bvec right_b_hat = "";
	
				if(u_hat.length() >= 2)
				{
					left_b_hat = odd(u_hat(0,length(u_hat)-2))+even(u_hat(0,length(u_hat)-2));
					right_b_hat = even(u_hat(0,length(u_hat)-2));
				}
				
				llr left = L_Nover2(L_1(0,L_1.length()/2-1),left_b_hat, Ls, done, (a-1)/2, c);
				llr right = L_Nover2(L_1(L_1.length()/2,L_1.length()-1),right_b_hat, Ls, done, (a-1)/2, c+L_1.length()/2);
				
				if(u_hat(u_hat.length()-1))
				{
	
					Ls(x,y) = right-left;
				}
				else
				{
					Ls(x,y) = right+left;
				}
			}
			done(x,y) = 1;
		}
			
		return Ls(x,y);
	}
}
