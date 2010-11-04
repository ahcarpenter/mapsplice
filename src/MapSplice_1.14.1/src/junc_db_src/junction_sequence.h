/*    
 *    junction_sequence.h		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Zeng Zheng
 *
 *    Authors: Zeng Zheng
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef JUNC_SEQ_H
#define JUNC_SEQ_H

#include <string>
using namespace std;

#pragma warning(disable:4996)

struct Junc_Seq
{
public:
	int start;
	string jump_code;
	string sequence;


	Junc_Seq()
	{}

	Junc_Seq(const Junc_Seq& my_junc_seq)
	{
		this->start=my_junc_seq.start;
		this->jump_code=my_junc_seq.jump_code;
		this->sequence=my_junc_seq.sequence;
	}

	~Junc_Seq()
	{}
};



#endif