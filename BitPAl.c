/**
*  Copyright (c) 2013,2019 Joshua Loving
*  Laboratory for Biocomputing and Informatics, Boston University
*  All rights reserved.
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*   + Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*   + Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*   + This source code may not be used in any program that is sold, any
*    derivative work must allow free distribution and modification.
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
*  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE LABORATORY FOR BIOCOMPUTING AND INFORMATICS
*  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


/*
*
*	This code implements a similarity distance using a scoring scheme
*	of match = 0, mismatch = -1, and gap = -1. This is analogous to the distance
*	score 0, 1, 1 of the Levenshtein distance. The values returned can be negated
* to obtain the Levenshtein distance.
*
* The algorithm implemented here for bit-parallel computation of similarity
* scores is described in https://doi.org/10.1093/bioinformatics/btu507
*
**/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define wordsize 64

int bitwise_alignment(char * s1, char* s2)
/**
*
* s1 and s2 are the strings to be aligned.
*
* The returned integer is the alignment score between s1 and s2.
*
*
**/
{

	int i, j, n, m;
	unsigned long long int bitmask;
	unsigned long long int Matches;
	unsigned long long int NotMatches;
	unsigned long long int all_ones = ~0x0000000000000000;
	unsigned long long int one = 0x0000000000000001;

	unsigned long long int DVpos1shift, DVzeroshift, DVneg1shift;
	unsigned long long int DHpos1, DHzero, DHneg1;
	unsigned long long int INITpos1s, INITzeros;
	unsigned long long int RemainDHneg1;
	unsigned long long int DHneg1tozero;
	unsigned long long int DVpos1shiftorMatch;
	unsigned long long int DVnot1to1shiftorMatch;
	unsigned long long int DHpos1orMatch;
	unsigned long long int add1;
	unsigned long long int add2;
	unsigned long long int add4;


	int score;


	char *iterate = s1;

  //Allocate memory for the matchvec, which encodes string1.
  //String1 is assumed to be ASCII encoded.
	unsigned long long int * matchvec = NULL;
	if (matchvec == NULL){
		matchvec = (unsigned long long int *) calloc(256, sizeof(unsigned long long int));
	}
	n = strlen (s1);
	m = strlen (s2);



	//*************************encode match strings A C G T N for string1
	//loop through string1 and store bits in matchA, matchC, etc.
	//position zero corresponds to column one in the score matrix, i.e., first character
	//so we start with i = 0 and bitmask = 1

	bitmask = 0x0000000000000001;
	for (i = 0; i < n; ++i)
	{
		matchvec[(*iterate)] |= bitmask;
		++iterate; bitmask <<= 1;
	}



	//intialize top row (penalty for initial gap, unless semi-global alignment is being done)
	DHneg1 = all_ones;
	DHzero = DHpos1 = 0;

	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		Matches = matchvec[*iterate];
		//Complement Matches
		NotMatches = ~Matches;

		//Finding the vertical values. 		//Find 1s
		INITpos1s = DHneg1 & Matches;
		DVpos1shift = (((INITpos1s + DHneg1) ^DHneg1) ^ INITpos1s);

		//set RemainingDHneg1
		RemainDHneg1 = DHneg1 ^ (DVpos1shift >> 1);
		//combine 1s and Matches
		DVpos1shiftorMatch = DVpos1shift | Matches;

		//set DVnot1to1shiftorMatch
		DVnot1to1shiftorMatch = ~(DVpos1shiftorMatch);


		//Find 0s
		DVzeroshift = ((DHzero & DVpos1shiftorMatch) | (DHneg1 & DVnot1to1shiftorMatch))<<1;
		//Find -1s
		DVneg1shift = all_ones ^ (DVpos1shift | DVzeroshift);
		//combine 1s and Matches
		DHpos1orMatch = DHpos1| Matches;
		//group -1topos0
		DHneg1tozero = all_ones^(DHpos1orMatch);
		//Find 0s
		DHzero = ((DVzeroshift & DHpos1orMatch) | (DVneg1shift & DHneg1tozero));
		//Find 1s
		DHpos1 = ((DVneg1shift & DHpos1orMatch) );
		//Find -1s
		DHneg1 = all_ones^(DHzero | DHpos1);
	}
	//find scores in last row
	add1 = DHzero;
	add2 = DHpos1;
	add4 = 0;

  //Reset bitmask to 1
	bitmask = 0x0000000000000001;
  //Initialize score to value in bottom left side of scoring matrix
	score = -1*m;

	for (i = 0; i < n; i++)
	{
		score += ((add1 & bitmask) >> i) * 1 + ((add2 & bitmask) >> i) * 2 + ((add4 & bitmask) >> i) * 4-1;
		//printf ("%4d", score);
		bitmask <<= 1;
	}


	free(matchvec);

	return score;
}




    int main(int argc, char * argv[]){
		char * fname;
		int i;
		if (argc == 1){

			fprintf(stderr, "usage: %s filename [-h]\n", argv[0]);
			exit(1);
		}
		for (i = 1; i < argc; i++)  // Skip program name
		{
			if (strcmp(argv[i], "-h") == 0)  /*if -h, print help */
			{
				printf(" \n");
				exit(0);
			}
			else
			{
				fname = argv[i];
			}
		}
        FILE * in = fopen(fname, "r");

        char *s1 = malloc(1000*sizeof(char));//One of the strings to align
        char *s2 = malloc(1000*sizeof(char));//The other string to align
        fgets(s1, 1000, in);
        while(fgets(s2, 1000, in) != NULL){


            s1[63] = '\0';
            s2[63] = '\0';
            printf ("%4d\t", bitwise_alignment(s1, s2, 1));
        }
        return 0;
        }
