/*    
 PolarCode provides functions for polar code training, encoding and decoder.
 Copyright (C) 2010  Robert G. Maunder
 
 This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
 The GNU General Public License can be seen at http://www.gnu.org/licenses/.
 */

/// \file main.cpp
/// \author Robert G Maunder
/// \date 06/08/2010
/// \version 1.1.0
/// \brief Trains polar codes and simulates their operation.


#include <cfloat>
#include "itpp/itcomm.h"
#include "PolarCode.h"

/// \ingroup polar
/// Trains polar codes and simulates their operation. 
/// Parameters are read from the command line and/or a configuration file. For example... 
/// \code./executable Training=1 N=1024 K=512 SNR=0 FrameCount=10000\endcode
/// This will train the information set \f$\mathcal{A}\f$ of an \f$(N,K)\f$ polar code for transmission over a BPSK-modulated AWGN channel having the specified SNR. The length and accuracy of the simulation depends on the value you use for FrameLength. The resultant parameters \f$(N,K,\mathcal{A})\f$ will be written into a configuration file having a filename of the form polar_N_K_SNR.
/// \code./executable Training=1 N=1024 ThresholdMI=0.999 SNR=0 FrameCount=10000\endcode 
/// This will train the information set \f$\mathcal{A}\f$ and automatically choose an appropriate information block length \f$K\f$ for a polar code having an encoded block length of \f$N\f$ used for transmission over a BPSK-modulated AWGN channel having the specified SNR. During this training process, the mutual information offered by each of the \f$N\f$ channels will be compared with ThresholdMI. If the mutual information is greater than ThresholdMI, then the channel will be admitted to the information set \f$\mathcal{A}\f$ and \f$K\f$ will be incremented. The length and accuracy of the simulation depends on the value you use for FrameLength. The resultant parameters \f$(N,K,\mathcal{A})\f$ will be written into a configuration file having a filename of the form polar_N_K_SNR.
/// \code./executable ConfigFile=polar_1024_512_+0.00 SNRStart=-10 SNRDelta=0.5 TargetBitCount=5120 TargetErrorCount=100\endcode
/// This will read the polar code parameters \f$(N,K,\mathcal{A})\f$ from the specified configuration file and simulate its transmission over a BPSK-modulated AWGN channel having the SNRs [SNRStart, SNRStart+SNRDelta, SNRStart+2SNRDelta, SNRStart+3SNRDelta, ...]. The simulation of each SNR will continue until at least TargetBitCount number of bits have been transmitted and TargetErrorCount number of errors have been observed. The results are written into a file having a filename of the form ConfigFile_SNRStart_ber. For each SNR, the results file provides the number of transmitted bits, the number of observed errors, the BER and the average number of Add, Compare and Select (ACS) operations performed per bit.
int main(int argc, char* argv[])
{
    // Create a parser and read parameters from the command line and the configuration file, if it exists
    Parser parser;
    parser.init(argc, argv);	
    
    // Should the polar code be trained or characterised?
    bool training = 0;
    if(parser.exist("Training"))
    {
        parser.get(training, "Training");	
    }
    
    // Train some parameters to optimise the polar code for a particular channel
    if(training)
    {
        // Allow the choice of whether to use the exact, lookup-table-aided or approximate Jacobian operator to be extracted from the parser. A default choice of the lookup-table-aided Jacobian operator is used if one cannot be extracted from the parser.
        unsigned int jacobian_mode = 1;
        if(parser.exist("Jacobian"))
        {
            parser.get(jacobian_mode, "Jacobian");
        }		
        if(jacobian_mode == 0)
        {
            jacobian_type = RobProb::Exact;
        }
        else if(jacobian_mode == 1)
        {
            jacobian_type = RobProb::Lookup;
        }
        else
        {
            jacobian_type = RobProb::Approx;
        }
        
        // Encoded block length
        int N = 8;
        if(parser.exist("N"))
        {
            parser.get(N, "N");	
        }
        
        // Specify the minimum mutual information that is regarded as 'good'
        double MI_threshold = 0.999;
        if(parser.exist("ThresholdMI"))
        {
            parser.get(MI_threshold, "ThresholdMI");	
        }
        
        // Information block length
        // Set K>0 to design a code having the coding rate K/N.
        // Set K=0 to design a code having a minimal coding rate that will give 'good' performance. The definition of 'good' depends on the choice of MI_threshold.
        int K = 0;
        if(parser.exist("K"))
        {
            parser.get(K, "K");	
        }
        
        // The more frames that are simulated, the better the results will be
        int frame_count = 10000;
        if(parser.exist("FrameCount"))
        {
            parser.get(frame_count, "FrameCount");	
        }
        
        // The trained parameters are dependent on the channel
        double snr = 0;
        if(parser.exist("SNR"))
        {
            parser.get(snr, "SNR");	
        }
        double N0 = 1.0/inv_dB(snr);
        
        // Use the histogram method to measure mutual information
        bool histogram = 1;
        if(parser.exist("Histogram"))
        {
            parser.get(histogram, "Histogram");	
        }
        
        
        // Use the histogram method to measure mutual information
        bool channel = 0;
        if(parser.exist("Rayleigh"))
        {
            parser.get(channel, "Rayleigh");	
        }
        
        bool bec_channel = 0;
        if(parser.exist("BEC"))
        {
            parser.get(bec_channel,"BEC");
        }
        
        bool bsc_channel = 0;
        if(parser.exist("BSC"))
        {
            parser.get(bsc_channel,"BSC");
        }
        
        double epsilon = 0.6;
        if(parser.exist("Epsilon"))
        {
            parser.get(epsilon, "Epsilon");
        }
        
        
        // Create the polar code, modulator and channel
        PolarCode polar_code(N);
        BPSK_c bpsk;
        AWGN_Channel awgn(N0);
        TDL_Channel rayleigh;
        
        // Storage for the simulation results
        bmat bits(frame_count, N);
        Mat<llr> llrs(frame_count, N);
        
        // Perform the simulation
        for(int frame_index = 0; frame_index < frame_count; frame_index++)
        {
            cout << frame_index << endl;
            // Randomise the frozen bits
            polar_code.randomise_u_A_c();
            //cout << "Frozen bits = " << polar_code.get_u_A_c() << endl;
            
            // Polar encoding
            bvec encoded_bits = polar_code.encode();
            //cout << "transmitted bits = " << encoded_bits << endl;
            
            llr_frame encoded_llrs(N);   
            
            //Choose the channel models - BSC,AWGN,BEC
            
            if(bsc_channel)
            {
                // BSC Channel 
                BSC bsc(epsilon);
                bvec bsc_bits = bsc(encoded_bits);
                //cout << "bsc bit = " <<  bsc_bits << endl; 
                
                //vec llr_double(N);
                for(int i = 0; i < N ; i++)
                {
                    if(bsc_bits(i) == 1)
                    {
                        if(epsilon >= 0.5)
                        {
                            encoded_llrs(i).set_to_plus_infinity();
                        }
                        else
                        {
                            encoded_llrs(i).set_to_minus_infinity();
                        }
                    }
                    else if(bsc_bits(i) == 0)
                    {
                        if(epsilon >= 0.5)
                        {
                            encoded_llrs(i).set_to_minus_infinity();
                        }
                        else
                        {
                            encoded_llrs(i).set_to_plus_infinity();
                        }
                    }   
                }
                //cout << "encoded_llrs of BSC = " << encoded_llrs << endl;
            }
            else if(bec_channel)
            {
                //BEC channel
                vec bec_bits(N);
                int erasure = 2;
                Uniform_RNG u(0.0,1.0);
                for(int i = 0; i < N; i++)
                {
                    if(u() <= epsilon) 
                    {
                        int temp_encoded = encoded_bits(i);
                        bec_bits(i) = temp_encoded + erasure;
                    }
                    else
                        bec_bits(i) = encoded_bits(i);
                }
                
               // cout << "received bits = " << bec_bits << endl;
                for(int i = 0; i < N ; i++)
                {
                    if(bec_bits(i) == 1)
                        encoded_llrs(i).set_to_minus_infinity();
                    else if(bec_bits(i) == 0)
                        encoded_llrs(i).set_to_plus_infinity();
                }
                //cout << "encoded_llrs of BEC = " << encoded_llrs << endl;   
            }
            else if (!bec_channel && !bsc_channel)
            {
                // BPSK modulation  
                // AWGN and AWGN+Rayleigh Channel
                cvec tx = bpsk.modulate_bits(encoded_bits);
                if(channel)
                {
                    Array<cvec> channel(1);
                    // Uncorrelated Rayleigh fading channel
                    cvec rx = awgn(rayleigh(tx,channel));
                    // Soft BPSK demodulation
                    encoded_llrs = to_llr_frame(bpsk.demodulate_soft_bits(rx,channel(0),N0));
                }
                else
                {
                    // AWGN channel
                    cvec rx = awgn(tx);
                    // Soft BPSK demodulation
                    encoded_llrs = to_llr_frame(bpsk.demodulate_soft_bits(rx,N0));
                }
                //cout<<"encoded_llrs of AWGN = " << encoded_llrs << endl;
            }
            
            // Polar decoding
            llr_frame decoded_llrs = polar_code.decode_frozen(encoded_llrs);
            cout<<"LLR frame of decoded bits" << decoded_llrs <<endl;
            
            // Store results
            bits.set_row(frame_index, polar_code.get_u_A_c()); //known frozen bits
            llrs.set_row(frame_index, decoded_llrs); //decoded frozen bits
        }
        
        // Measure mutual informations of each channel
        vec MIs(N);
        for(int i = 0; i < N; i++)
        {
            if(histogram)
            {
                MIs(i) = mutual_information(llrs.get_col(i), bits.get_col(i));
            }
            else
            {
                MIs(i) = mutual_information(llrs.get_col(i));
            }
        }
        
        // Determine and display the best channels
        ivec A = sort_index(MIs)+1;
        sort(MIs);
        cout << "Channel MI" << endl;
        cout << MIs << endl;
        //cout << "Channel index " <<endl;
        //cout << A << endl;
        //cout << "MI threshold " << MI_threshold <<endl;
        
        if(K==0) //this is only done if K is not specified. Otherwise, the number of good channels = K
        {
            // See which channels are 'good'
            while(MIs(K) < MI_threshold)
            {
                K++;
            }
        }
        
        // Select the best channels
        A=A(N-K,N-1);
        sort(A);
        
        // Output parameters which can be used in a BER simulation
        cout << "N = " << N << endl;
        cout << "K = " << K << endl;
        cout << "A = " << A << endl;
        cout << "No. of channels = " << A.length() << endl;
        cout << endl << endl << endl;
        
        
        // Choose a file to save the results into.
        char output_filename[100];
        sprintf(output_filename, "polar_%d_%d_%+0.2f", N, K, snr);
        fstream fout;
        fout.open(output_filename, ios::out | ios::trunc);
        /*fout << "Channel MI = " << endl;
         for(int k=0;k<MIs.length()-1;k++)
         {
         fout << MIs(k) << endl;
         }*/
        
        
        fout << "N=" << N << endl;
        fout << "K=" << K << endl;
        fout << "A=" << '[';
        for(int k = 0; k<A.length()-1; k++)
        {
            fout << A(k) << ',';
        }
        fout << A(A.length()-1) << ']' << endl;	
        
	fout << "ChannelModel=";
        if(bsc_channel)
            fout << "1" <<endl;
        else if(bec_channel)
            fout << "2" <<endl;
        else if (!bsc_channel &&  !bec_channel)
            fout<< "3" <<endl;
        
        fout << "Epsilon=" << epsilon << endl;
        fout.close();
        
        
        
    }
    else // Characterise the polar code by drawing its BER plot
    {
        string config_file = "polar";	
        if(parser.exist("ConfigFile"))
        {
            config_file = parser.get_string("ConfigFile");	
        }
        fstream fin;
        fin.open(config_file.c_str(), ios::in);
        if(!fin.fail())
        {
            parser.init(config_file, argc, argv);
            fin.close();
        }
	
        // Allow the choice of whether to use the exact, lookup-table-aided or approximate Jacobian operator to be extracted from the parser. A default choice of the lookup-table-aided Jacobian operator is used if one cannot be extracted from the parser.
        unsigned int jacobian_mode = 1;
        if(parser.exist("Jacobian"))
        {
            parser.get(jacobian_mode, "Jacobian");
        }		
        if(jacobian_mode == 0)
        {
            jacobian_type = RobProb::Exact;
        }
        else if(jacobian_mode == 1)
        {
            jacobian_type = RobProb::Lookup;
        }
        else
        {
            jacobian_type = RobProb::Approx;
        }
        
        // Encoded block length
        int N=8;
        if(parser.exist("N"))
        {
            parser.get(N, "N");	
        }
        
        // Information block length
        int K=5;
        if(parser.exist("K"))
        {
            parser.get(K, "K");	
        }
        
        // Information set
        ivec A="3,5,6,7,8";		
        if(parser.exist("A"))
        {
            parser.get(A, "A");	
        }
        
        int channelModel = 3;
        if(parser.exist("ChannelModel"))
        {
            parser.get(channelModel,"ChannelModel");
        }
        
        double epsilon = 0.6;
        if(parser.exist("Epsilon"))
        {
            parser.get(epsilon,"Epsilon");
        }
        
        // Set the default starting Signal to Noise Ratio (SNR).
        double snr_start = -10.0; // dB
        // Read in the starting SNR.
        if(parser.exist("SNRStart"))
        {
            parser.get(snr_start, "SNRStart");	
        }
	
        // Set the default SNR delta.
        double snr_delta = 0.5; // dB
        // Read in the SNR delta.
        if(parser.exist("SNRDelta"))
        {
            parser.get(snr_delta, "SNRDelta");	
        }
	
        // Set the default stopping SNR.
        double snr_stop = DBL_MAX; // dB
        // Read in the starting SNR.
        if(parser.exist("SNRStop"))
        {
            parser.get(snr_stop, "SNRStop");	
        }
        
        // Set the default stopping BER.
        double ber_stop = 0;
        // Read in the starting SNR.
        if(parser.exist("BERStop"))
        {
            parser.get(ber_stop, "BERStop");	
        }
	
        // Minimum number of bits to simulate at each SNR
        int target_bit_count = K*100;
        if(parser.exist("TargetBitCount"))
        {
            parser.get(target_bit_count, "TargetBitCount");	
        }
        
        // Minimum number of errors to observe at each SNR
        int target_error_count = 100;
        if(parser.exist("TargetErrorCount"))
        {
            parser.get(target_error_count, "TargetErrorCount");	
        }
        
        // Use the histogram method to measure mutual information
        bool channel = 1;
        if(parser.exist("Rayleigh"))
        {
            parser.get(channel, "Rayleigh");	
        }
        
        // Choose a file to save the results into.
        char output_filename[100];
        sprintf(output_filename, "%s_%+0.2f_ber", config_file.c_str(), snr_start);
	
        // Open the output file.
        fstream fout;
        fout.open(output_filename, ios::out | ios::trunc);
        
        // Output the header information.
        fout << "!SNR\tberrors\tbits\tBER\tferrors\tframes\tFER\tACS" << endl;
        cout << "SNR\tberrors\tbits\tBER\tferrors\tframes\tFER\tACS" << endl;
        
        
        // Create the polar code, modulator and channel
        PolarCode polar_code(N,K,A);
        BPSK_c bpsk;
        AWGN_Channel awgn;
        TDL_Channel rayleigh;
        BSC bsc(epsilon);
        int erasure = 2;
        Uniform_RNG u(0.0,1.0);
        
        // Loop until the job is killed or until the SNR or BER target is reached.
        double ber = 1.0;
        double fer = 1.0;
        for(double snr = snr_start; snr < snr_stop && ber > ber_stop; snr += snr_delta)
        {
            // Convert from SNR (in dB) to noise power spectral density.
            double N0 = 1.0/inv_dB(snr);
            awgn.set_noise(N0);
            
            int error_count = 0;
            int bit_count = 0;
            int total_acs_operations = 0;
            int frame_error_count = 0;
            int frame_count = 0;
            
            // Keep going until enough errors have been observed. This runs the simulation only as long as is required to keep the BER vs SNR curve smooth.
            while(bit_count < target_bit_count || error_count < target_error_count)
            {
                // Generate some information bits;
                bvec information_bits = randb(K);
                
                // Randomise the frozen bits
                polar_code.randomise_u_A_c();
                
                // Polar encoding
                bvec encoded_bits = polar_code.encode(information_bits);
                //cout << "information bits = " << information_bits << endl;
                //cout << "encoded transmitted bits = " << encoded_bits <<endl;
                
                // BPSK modulation
                
                llr_frame encoded_llrs(N);
                
                //choose channel
                // BSC Channel 
                
                if(channelModel == 1)
                {
                    bvec bsc_bits = bsc(encoded_bits);
                    
                    // BEC channel LLR decoding
                    for(int i = 0; i < N ; i++)
                    {
                        if(bsc_bits(i) == 1)
                        {
                            if(epsilon >= 0.5)
                            {
                                encoded_llrs(i).set_to_plus_infinity();
                            }
                            else
                            {
                                encoded_llrs(i).set_to_minus_infinity();
                            }
                        }
                        else if(bsc_bits(i) == 0)
                        {
                            if(epsilon >= 0.5)
                            {
                                encoded_llrs(i).set_to_minus_infinity();
                            }
                            else
                            {
                                encoded_llrs(i).set_to_plus_infinity();
                            }
                        }   
                    }
                }
                
                //BEC channel
                else if(channelModel == 2)
                {
                    // Channel Modelling
                    vec bec_bits(N);
                    
                    for(int i = 0; i < N; i++)
                    {
                        if(u() <= epsilon) 
                        {
                            int temp_encoded = encoded_bits(i);
                            bec_bits(i) = temp_encoded + erasure;
                        }
                        else
                            bec_bits(i) = encoded_bits(i);
                    }
                    //cout << "bec bits" << bec_bits << endl;
                    //getchar();
                    
                    // BEC LLR decoding
                    for(int i = 0; i < N ; i++)
                    {
                        if(bec_bits(i) == 1)
                            encoded_llrs(i).set_to_minus_infinity();
                        else if(bec_bits(i) == 0)
                            encoded_llrs(i).set_to_plus_infinity();
                    }   
                }
                
                // AWGN channel
                if(channelModel == 3)
                {
                    cvec tx = bpsk.modulate_bits(encoded_bits);
                    if(channel)
                    {
                        Array<cvec> channel(1);
                        // Uncorrelated Rayleigh fading channel
                        cvec rx = awgn(rayleigh(tx,channel));
                        // Soft BPSK demodulation
                        encoded_llrs = to_llr_frame(bpsk.demodulate_soft_bits(rx,channel(0),N0));
                    }
                    else
                    {
                        // AWGN channel
                        cvec rx = awgn(tx);
                        // Soft BPSK demodulation
                        encoded_llrs = to_llr_frame(bpsk.demodulate_soft_bits(rx,N0));
                    }
                }
                
                // Polar decoding
                int acs_operations;
                llr_frame decoded_llrs = polar_code.decode(encoded_llrs, acs_operations);
                //cout << "encoded received bits = " << hard(encoded_llrs) << endl;
                //cout << "decoded bits = " << hard(decoded_llrs) << endl;
                // Update the results
                error_count += sum(to_ivec(information_bits + hard(decoded_llrs)));
                //cout << "Error count = " << error_count << endl;
                if(sum(to_ivec(information_bits + hard(decoded_llrs))) > 0)
                {
                    frame_error_count++;
                }
                bit_count += A.length();
                frame_count++;
                total_acs_operations += acs_operations;
                
                //cout<<endl<<endl;
                //getchar(); //debug point
            }
            
            // Calculate the BER.
            ber = static_cast<double>(error_count)/static_cast<double>(bit_count);
            fer = static_cast<double>(frame_error_count)/static_cast<double>(frame_count);
            
            // Output the results.
            fout << snr << '\t' << error_count << '\t' << bit_count << '\t' << ber << '\t' << frame_error_count << '\t' << frame_count << '\t' << fer << '\t' << static_cast<double>(total_acs_operations)/static_cast<double>(bit_count) << endl;
            cout << snr << '\t' << error_count << '\t' << bit_count << '\t' << ber << '\t' << frame_error_count << '\t' << frame_count << '\t' << fer << '\t' << static_cast<double>(total_acs_operations)/static_cast<double>(bit_count) << endl;
        }
        
        // Close the output file.
        fout.close();
        
    }
    return 0;
}

