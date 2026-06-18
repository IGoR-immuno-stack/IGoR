/*
 * main.cpp - IGoR Demo standalone application
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  This standalone application extracts the demo code from the main IGoR application.
 *  It demonstrates:
 *  - Reading TCRb genomic templates
 *  - Aligning sequences to templates
 *  - Creating a TCRb model with error rate and marginals
 *  - Reading and writing model parameters
 *  - Inferring a model from sequences (10 iterations of EM)
 *  - Generating sequences from the obtained model
 */

#include <igor/Core/Deletion.h>
#include <igor/Core/Insertion.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Model_Parms.h>
#include <igor/Core/Rec_Event.h>
#include <igor/Core/Singleerrorrate.h>
#include <igor/Core/Model_marginals.h>
#include <igor/Core/Aligner.h>
#include <igor/Core/GenModel.h>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/Counter.h>
#include <igor/Core/Coverageerrcounter.h>
#include <igor/Core/Bestscenarioscounter.h>
#include <igor/Core/Pgencounter.h>
#include <igor/Core/Errorscounter.h>
#include <igor/Core/Utils.h>

#include <iostream>
#include <chrono>
#include <string>

using namespace std;

std::string IGOR_DATA_DIR = "../..";

int main(int argc, char* argv[]) {

	/*Run this sample demo code
	 *
	 * Outline:
	 *
	 * Read TCRb genomic templates
	 *
	 * Align the sequences contained in the /demo/murugan_naive1_noncoding_demo_seqs.txt file to those templates
	 *
	 * Create a TCRb model, a simple error rate and the corresponding marginals
	 *
	 * Show reads and write functions for the model and marginals
	 *
	 * Infer a model from the sequences (perform 10 iterations of EM)
	 *
	 * Generate sequences from the obtained model
	 *
	 */

	// Working directory path (can be passed as command line argument)
	string cl_path;
	if (argc > 1) {
		cl_path = string(argv[1]);
		if (cl_path[cl_path.size()-1] != '/') {
			cl_path += "/";
		}
	} else {
		cl_path = "./";
	}

	clog<<"Reading genomic templates"<<endl;

	vector<pair<string,string>> v_genomic = read_genomic_fasta( string(IGOR_DATA_DIR) + "/demo/genomicVs_with_primers.fasta");

	vector<pair<string,string>> d_genomic = read_genomic_fasta( string(IGOR_DATA_DIR) + "/demo/genomicDs.fasta");

	vector<pair<string,string>> j_genomic = read_genomic_fasta( string(IGOR_DATA_DIR) + "/demo/genomicJs_all_curated.fasta");

	//Declare substitution matrix used for alignments(nuc44 here)
	double nuc44_vect [] = { // A,C,G,T,R,Y,K,M,S,W,B,D,H,V,N
		5,-14,-14,-14,-14,2,-14,2,2,-14,-14,1,1,1,0,
		-14,5,-14,-14,-14,2,2,-14,-14,2,1,-14,1,1,0,
		-14,-14,5,-14,2,-14,2,-14,2,-14,1,1,-14,1,0,
		-14,-14,-14,5,2,-14,-14,2,-14,2,1,1,1,-14,0,
		-14,-14,2,2,1.5,-14,-12,-12,-12,-12,1,1,-13,-13,0,
		2,2,-14,-14,-14,1.5,-12,-12,-12,-12,-13,-13,1,1,0,
		-14,2,2,-14,-12,-12,1.5,-14,-12,-12,1,-13,-13,1,0,
		2,-14,-14,2,-12,-12,-14,1.5,-12,-12,-13,1,1,-13,0,
		2,-14,2,-14,-12,-12,-12,-12,1.5,-14,-13,1,-13,1,0,
		-14,2,-14,2,-12,-12,-12,-12,-14,1.5,1,-13,1,-13,0,
		-14,1,1,1,1,-13,1,-13,-13,1,0.5,-12,-12,-12,0,
		1,-14,1,1,1,-13,-13,1,1,-13,-12,0.5,-12,-12,0,
		1,1,-14,1,-13,1,-13,1,-13,1,-12,-12,0.5,-12,0,
		1,1,1,-14,-13,1,1,-13,1,-13,-12,-12,-12,0.5,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Matrix<double> nuc44_sub_matrix(15,15,nuc44_vect);

	//Instantiate aligner with substitution matrix declared above and gap penalty of 50
	Aligner v_aligner = Aligner(nuc44_sub_matrix , 50 , V_gene);
	v_aligner.set_genomic_sequences(v_genomic);

	Aligner d_aligner = Aligner(nuc44_sub_matrix , 50 , D_gene);
	d_aligner.set_genomic_sequences(d_genomic);

	Aligner j_aligner (nuc44_sub_matrix , 50 , J_gene);
	j_aligner.set_genomic_sequences(j_genomic);

	clog<<"Reading sequences and aligning"<<endl;
	typedef std::chrono::system_clock myclock;
	myclock::time_point begin_time, end_time;

	begin_time = myclock::now();


	vector<pair<const int, const string>> indexed_seqlist = read_txt( string(IGOR_DATA_DIR) + "/demo/murugan_naive1_noncoding_demo_seqs.txt" ); //Could also read a FASTA file <code>read_fasta()<\code> or indexed sequences <code>read_indexed_seq_csv()<\code>

	cl_path+="igor_demo/";
	system(&("mkdir " + cl_path )[0]);

	v_aligner.align_seqs( string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_V.csv"),indexed_seqlist,50,true,INT16_MIN,-155);
	//v_aligner.write_alignments_seq_csv(path + string("alignments_V.csv") , v_alignments);

	d_aligner.align_seqs(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_D.csv"),indexed_seqlist,0,false);
	//d_aligner.write_alignments_seq_csv(path + string("alignments_D.csv") , d_alignments);

	j_aligner.align_seqs(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_J.csv"),indexed_seqlist,10,true,42,48);

	end_time= myclock::now();
	chrono::duration<double> elapsed = end_time - begin_time;
	clog<<"Alignments procedure lasted: "<<elapsed.count()<<" seconds"<<endl;
	clog<<"for "<<indexed_seqlist.size()<<" TCRb sequences of 60bp(from murugan and al), against ";
	clog<<v_genomic.size()<<" Vs,"<<d_genomic.size()<<" Ds, and "<<j_genomic.size()<<" Js full sequences"<<endl;

	//unordered_map<int,forward_list<Alignment_data>> j_alignments = j_aligner.align_seqs(indexed_seqlist,10,true,42,48);
	//j_aligner.write_alignments_seq_csv(path + string("alignments_J.csv") , j_alignments);


	write_indexed_seq_csv(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_indexed_seq.csv") , indexed_seqlist);



	clog<<"Construct the model"<<endl;
	//Construct a TCRb model
	Gene_choice v_choice(V_gene,v_genomic);
	v_choice.set_nickname("v_choice");
	v_choice.set_priority(7);
	Gene_choice d_choice(D_gene,d_genomic);
	d_choice.set_nickname("d_gene");
	d_choice.set_priority(6);
	Gene_choice j_choice(J_gene,j_genomic);
	j_choice.set_nickname("j_choice");
	j_choice.set_priority(7);

	Deletion v_3_del(V_gene_seq,Three_prime,make_pair(-4,16));//16
	v_3_del.set_nickname("v_3_del");
	v_3_del.set_priority(5);
	Deletion d_5_del(D_gene_seq,Five_prime,make_pair(-4,16));
	d_5_del.set_nickname("d_5_del");
	d_5_del.set_priority(5);
	Deletion d_3_del(D_gene_seq,Three_prime,make_pair(-4,16));
	d_3_del.set_nickname("d_3_del");
	d_3_del.set_priority(5);
	Deletion j_5_del(J_gene_seq,Five_prime,make_pair(-4,18));
	j_5_del.set_nickname("j_5_del");
	j_5_del.set_priority(5);

	Insertion vd_ins(VD_ins_seq,make_pair(0,30));
	vd_ins.set_nickname("vd_ins");
	vd_ins.set_priority(4);
	Insertion dj_ins(DJ_ins_seq,make_pair(0,30));
	dj_ins.set_nickname("dj_ins");
	dj_ins.set_priority(2);

	Dinucl_markov markov_model_vd(VD_ins_seq);
	markov_model_vd.set_nickname("vd_dinucl");
	markov_model_vd.set_priority(3);

	Dinucl_markov markov_model_dj(DJ_ins_seq);
	markov_model_dj.set_nickname("dj_dinucl");
	markov_model_dj.set_priority(1);


	Model_Parms parms;

	//Add nodes to the graph
	parms.add_event(&v_choice);
	parms.add_event(&d_choice);
	parms.add_event(&j_choice);

	parms.add_event(&v_3_del);
	parms.add_event(&d_3_del);
	parms.add_event(&d_5_del);
	parms.add_event(&j_5_del);

	parms.add_event(&vd_ins);
	parms.add_event(&dj_ins);

	parms.add_event(&markov_model_vd);
	parms.add_event(&markov_model_dj);


	//Add correlations
	parms.add_edge(&v_choice,&v_3_del);
	parms.add_edge(&j_choice,&j_5_del);
	parms.add_edge(&d_choice,&d_3_del);
	parms.add_edge(&d_choice,&d_5_del);
	parms.add_edge(&d_5_del,&d_3_del);
	parms.add_edge(&j_choice,&d_choice);


	//Create the corresponding marginals
	Model_marginals model_marginals(parms);
	model_marginals.uniform_initialize(parms); //Can also start with a random prior using random_initialize()

	//Instantiate an error rate
	Single_error_rate error_rate(0.001);

	parms.set_error_ratep(&error_rate);

	clog<<"Write and read back the model"<<endl;
	//Write the model_parms into a file
	parms.write_model_parms(string(cl_path + "/demo_write_model_parms.txt"));

	//Write the marginals into a file
	model_marginals.write2txt(string(cl_path + "/demo_write_model_marginals.txt"),parms);

	//Read a model and marginal pair
	Model_Parms read_model_parms;
	read_model_parms.read_model_parms(string(cl_path + "/demo_write_model_parms.txt"));
	Model_marginals read_model_marginals(read_model_parms);
	read_model_marginals.txt2marginals(string(cl_path + "/demo_write_model_marginals.txt"),read_model_parms);

	//Instantiate a Counter
	map<size_t,shared_ptr<Counter>> counters_list;
	 //Collect gene coverage and errors
	shared_ptr<Counter> coverage_counter_ptr(new Coverage_err_counter(cl_path + "/run_demo/",V_gene,1,false,false));
	counters_list.emplace(0,coverage_counter_ptr);

	 //Collect 10 best scenarios per sequence during the last iteration
	shared_ptr<Counter>best_sc_ptr(new Best_scenarios_counter(10 , cl_path + "/run_demo/" ,true ));
	counters_list.emplace(1,best_sc_ptr);

	 //Collect sequence generation probability during last iteration
	shared_ptr<Counter> pgen_counter_ptr(new Pgen_counter (cl_path + "/run_demo/"));
	counters_list.emplace(2,pgen_counter_ptr);

	shared_ptr<Counter> errors_counter(new Errors_counter (10,string(cl_path + "/run_demo/")));
	counters_list.emplace(3,errors_counter);

	//Instantiate the high level GenModel class
	//This class allows to make most useful high level operations(model inference/Pgen computation , sequence generation)
	GenModel gen_model(read_model_parms,read_model_marginals,counters_list);


	//Inferring a model

	//Read alignments
	//vector<pair<const int, const string>> indexed_seqlist = read_indexed_csv(path+ string(argv[2]) + string("indexed_seq.csv"));
	unordered_map<int,pair<string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments = read_alignments_seq_csv_score_range(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_V.csv"), V_gene , 55 , false , indexed_seqlist  );//40//35
	sorted_alignments = read_alignments_seq_csv_score_range(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_D.csv"), D_gene , 35 , false , indexed_seqlist , sorted_alignments);//30//15
	sorted_alignments = read_alignments_seq_csv_score_range(string(cl_path + "/murugan_naive1_noncoding_demo_seqs") + string("_alignments_J.csv"), J_gene , 10 , false , indexed_seqlist , sorted_alignments);//30//20

	vector<tuple<int,string,unordered_map<Gene_class,vector<Alignment_data>>>> sorted_alignments_vec = map2vect(sorted_alignments);

	//Infer the model
	clog<<"Infer model"<<endl;

	begin_time = myclock::now();
	system(&("mkdir " + cl_path + "run_demo")[0]);
	gen_model.infer_model(sorted_alignments_vec , 4 , string(cl_path + "/run_demo/") , true ,1e-35,0.0001);

	end_time= myclock::now();
	elapsed = end_time - begin_time;
	clog<<"Model inference procedure lasted: "<<elapsed.count()<<" seconds"<<endl;


	//Generate sequences
	clog<<"Generate sequences"<<endl;
	auto generated_seq =  gen_model.generate_sequences(100,false);//Without errors
	gen_model.write_seq_real2txt(string(cl_path + "/generated_seqs_indexed_demo.csv"), string(cl_path + "/generated_seqs_realizations_demo.csv") , generated_seq);//Member function will be changed

	clog<<"Demo completed successfully!"<<endl;
	return EXIT_SUCCESS;
}
