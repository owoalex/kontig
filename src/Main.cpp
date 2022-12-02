/*

Copyright (c) Alex Baldwin 2022.

This file is part of kontig.

kontig is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

kontig is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with kontig. If not, see <https://www.gnu.org/licenses/>. 

*/

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

#include <stdlib.h>
#include <cstring>

#include <stdio.h>

#include "Read.h"
#include "KMer.h"
#include "KMerTreeBuilder.h"
#include "ProposedContig.h"
#include "Contig.h"

int main(int argc, char** argv) {
    std::cout << "kontig DNA read assembler\n\n";
    char* input_file = NULL;
    char* output_file = NULL;
    uint8_t kontig_mode = 0;
    int kmer_length = 32; // 4^n = 2^2n; 128 length = 256 bits of information per kmer;
    /*
     * 0: Display help and exit
     * 1: Generate map from FASTQ file
     * 2: 
     * 
     * 
     * 
     */
    
    for (int i = 1; i < argc; ++i) {
        char* arg = argv[i];
        if (strcmp(arg,"--input") == 0 || strcmp(arg,"-i") == 0) {
            i++;
            input_file = argv[i];
            continue;
        }
        if (strcmp(arg,"--output") == 0 || strcmp(arg,"-o") == 0) {
            i++;
            output_file = argv[i];
            continue;
        }
        if (strcmp(arg,"--ksize") == 0) {
            i++;
            kmer_length = atoi(argv[i]);
            continue;
        }
        if (strcmp(arg,"help") == 0 || strcmp(arg,"-h") == 0 || strcmp(arg,"--help") == 0) {
            kontig_mode = 0;
            continue;
        }
        if (strcmp(arg,"genmap") == 0) {
            kontig_mode = 1;
            continue;
        }
        if (strcmp(arg,"pathfind") == 0) {
            kontig_mode = 2;
            continue;
        }
        std::cout << "Unknown command " << arg << "\n\n";
        kontig_mode = 0;
        break;
        //bootVmImmediate = arg;
    }

    if (kontig_mode == 0) {
        std::cout << "General command structure\n";
        std::cout << "kontig <mode> [options and flags]\n";
        std::cout << "\n";
        std::cout << "MODES:\n";
        std::cout << "\n";
        std::cout << "genmap - Generate .krm (kontig read map) of kmers derives from FASTQ files and connect them.\n";
        std::cout << "[TODO] concat - Consolidates multiple .krm files into one, useful when batch processing on distributed systems.\n";
        std::cout << "[TODO] pathfind - Find valid paths through the map to generate .krf (kontig route fragment) files.\n";
        std::cout << "[TODO] combine - Consolidates multiple .krf files into one, useful when batch processing on distributed systems.\n";
        std::cout << "[TODO] export - Turns .krf (kontig route fragment) files into FASTA files for further processing by other tools.\n";
        std::cout << "\n";
        std::cout << "OPTIONS:\n";
        std::cout << "\n";
        std::cout << "--input <file>\n";
        std::cout << "-i <file>\n";
        std::cout << "Input file(s) for the given mode. If specifying multiple files, this flag must be used multiple times.\n";
        std::cout << "\n";
        std::cout << "--output <file>\n";
        std::cout << "-o <file>\n";
        std::cout << "Output file for the given mode.\n";
        std::cout << "\n";
        std::cout << "--ksize <integer>\n";
        std::cout << "Set the K-Mer size used for map generation in genmap mode. 32 is the default value if not specified.\n";
        std::cout << "\n";
        exit(0);
    }
    
    //if (input_file == NULL) {
    //    std::cout << "Need to specify input file with '-i'.\n";
    //    exit(1);
    //}

    
    if (kontig_mode == 1) {
        
        uint64_t prog = 0;
        
        std::ifstream input_stream(input_file);
        
        std::ofstream output_stream(output_file);
        
        char* input_buffer = (char*) malloc(1024);

        input_stream.seekg(0, input_stream.end);
        int input_length = input_stream.tellg();
        input_stream.seekg(0, input_stream.beg);
        
        if (input_length == -1) {
            printf("Could not read input file.\n");
            exit(1);
        }
        
        printf("About to read %d bytes\n", input_length);
        
        bool still_readable = true;
        input_stream.read(input_buffer, 1024);
        int chars_read = input_stream.gcount();
        
        char* total_input = (char*) malloc(input_length);
        
        int input_position = 0;
        do {
            //*output_stream << "test";
            for (int i = 0; i < chars_read; i++) {
                if (input_buffer[i] == '\r') {
                    input_length--;
                } else {
                    total_input[input_position] = input_buffer[i];
                    input_position++;
                } 
            }
            
            prog++;
            if (prog > (1024 * 8)) {
                printf("\33[2K\r");
                printf("%.3f", (((double) input_position) / ((double) input_length)) * 100.0);
                std::cout << "%" << std::flush;
                prog = 0;
            }
            
            input_stream.read(input_buffer, 1024);
            chars_read = input_stream.gcount();
            still_readable = chars_read > 0;
        } while (still_readable);
        
        free(input_buffer);
        
        // Now parse the file
        
        input_position = 0;
        
        char* read_buffer = (char*) malloc(1024 * 64);
        char* name;
        char* sequence;
        char* quality;
        int read_buffer_position = 0;
        int mode = 0;
        std::vector<Read*> reads;
        /*
        * 0 : Scanning for separator
        * 1 : Reading in name
        * 2 : Reading in genome
        * 3 : Reading in quality
        */
        
        while (input_position < input_length) {
            switch (mode) {
                case 3:
                    if (total_input[input_position] == '\n') {
                    } else if (total_input[input_position] == '@') {
                        mode = 1;
                        quality = (char*) malloc(read_buffer_position + 2); // Make a new string of length of the string
                        std::memcpy(quality,read_buffer,read_buffer_position);
                        quality[read_buffer_position] = 0;
                        read_buffer_position = 0;
                        input_position++;
                        reads.push_back(new Read(name, sequence, quality));
                    } else {
                        read_buffer[read_buffer_position] = total_input[input_position];
                        read_buffer_position++;
                    }
                    break;
                case 2:
                    if (total_input[input_position] == '\n') {
                    } else if (total_input[input_position] == '+') {
                        mode = 3;
                        sequence = (char*) malloc(read_buffer_position + 2); // Make a new string of length of the string
                        std::memcpy(sequence,read_buffer,read_buffer_position);
                        sequence[read_buffer_position] = 0;
                        read_buffer_position = 0;
                    } else {
                        read_buffer[read_buffer_position] = total_input[input_position];
                        read_buffer_position++;
                    }
                    break;
                case 1:
                    if (total_input[input_position] == '\n') {
                        mode = 2; // Start reading genome
                        name = (char*) malloc(read_buffer_position + 2); // Make a new string of length of the string
                        std::memcpy(name,read_buffer,read_buffer_position);
                        name[read_buffer_position] = 0;
                        read_buffer_position = 0;
                    } else {
                        read_buffer[read_buffer_position] = total_input[input_position];
                        read_buffer_position++;
                    }
                    break;
                case 0:
                default:
                    if (total_input[input_position] == '@') {
                        mode = 1;
                    }
                    break;
            }
            input_position++;
        }
        
        free(total_input);
        printf("\rRead %d characters from FASTQ file\n", input_position);
        
        // process below
        
        KMerSet** kmer_set_set = (KMerSet**) malloc(sizeof(KMerSet*) * reads.size());
        int total_kmers = 0;
        prog = 0;
        
        printf("About to generate K-Mers of length %d\n", kmer_length);
        
        for (std::vector<Read*>::size_type i = 0; i < reads.size(); i++) {
            struct KMerSet* kms = reads[i]->generateKmers(kmer_length);
            total_kmers += kms->length;
            kmer_set_set[i] = kms;
            prog++;
            if (prog > (1024 * 8)) {
                printf("\33[2K\r");
                printf("%.3f", (((double) i) / ((double) reads.size())) * 100.0);
                std::cout << "%" << std::flush;
                prog = 0;
            }
        }
        
        printf("\rGenerated %d K-Mers\n", total_kmers);
        printf("About to sort K-Mers into tree\n");
        
        //KMer** kmers = (KMer**) malloc(sizeof(KMer*) * total_kmers);
        
        //int current_kmer = 0;
        
        
        KMerTreeBuilder* kmer_tree = new KMerTreeBuilder();
        prog = 0;
        
        for (std::vector<Read*>::size_type i = 0; i < reads.size(); i++) {
            for (int j = 0; j < kmer_set_set[i]->length; j++) {
                kmer_tree->insertKmer(kmer_set_set[i]->kmers[j]);
                
                prog++;
                if (prog > (1024 * 256)) {
                    printf("\33[2K\r");
                    printf("%.3f", (((double) i) / ((double) reads.size())) * 100.0);
                    std::cout << "%" << std::flush;
                    prog = 0;
                }
            }
        }
        
        //free(kmer_set_set);
        printf("\rK-Mers loaded into tree\n");
        
        // Now connect KMers

        printf("Connecting K-Mers\n");
        
        uint64_t total_connections = 0;
        uint16_t minimum_distance = 32;
        uint16_t maximum_mismatch = 2;
        prog = 0;
        std::vector<KMerEdge*> kmer_graph;
        
        for (std::vector<Read*>::size_type i = 0; i < reads.size(); i++) {
            KMer* head_kmer = reads[i]->kmers->kmers[0];
            //KMer* tail_kmer = reads[i]->kmers->kmers[reads[i]->kmers->length];
            //head_kmer->quickRef;
            KMerSortingTreeNode* cnode = kmer_tree->root;
            for (int j = 0; j < 64; j++) {
                if ((head_kmer->quickref >> j) & 0b1) {
                    cnode = cnode->high;
                } else {
                    cnode = cnode->low;
                }
            }
            KMerLinkedList* ll = cnode->values;
            while (ll != nullptr) {
                KMer* candidate_kmer = ll->value;
                if (head_kmer->source != candidate_kmer->source) {
                    if (candidate_kmer->offset > minimum_distance) { // Let's not bother with kmers that won't extend much
                        if (head_kmer->getMismatchedBases(candidate_kmer) < maximum_mismatch) {
                            int64_t match_quality = head_kmer->getMatchQuality(candidate_kmer);
                            
                            KMerEdge* kmer_edge = new KMerEdge();
                            kmer_edge->src = candidate_kmer;
                            kmer_edge->ext = head_kmer;
                            kmer_edge->weight = match_quality;
                            kmer_graph.push_back(kmer_edge);
                            //candidate_kmer->kmerEdgesForward.push_back(kmer_edge);
                            //head_kmer->kmerEdgesBackward.push_back(kmer_edge);
                            candidate_kmer->source->kmerEdgesForward.push_back(kmer_edge);
                            head_kmer->source->kmerEdgesBackward.push_back(kmer_edge);
                            
                            candidate_kmer->source->usedNTimes++;
                            
                            total_connections++;
                            prog++;
                            if (prog > (1024 * 64)) {
                                printf("\33[2K\r");
                                printf("%.3f", (((double) i) / ((double) reads.size())) * 100.0);
                                std::cout << "%" << std::flush;
                                //printf("%s EXTENDS TO\n%s\n", candidate_kmer->sequence, head_kmer->sequence);
                                prog = 0;
                            }
                        }
                    }
                }
                ll = ll->next;
            }
        }
        
        printf("\rMade %ld connections\n", total_connections);
        
        // Connections made, time to make proposed contigs via greedy algorithm
        
        prog = 0;
        uint64_t total_contigs = 0;
        
        //uint64_t starting_edge = 32; // This will be random eventually for multithreading
        KMerEdge* current_edge;
        ProposedContig* current_contig;
        Contig* generated_contig;
        
        for (std::vector<KMerEdge*>::size_type i = 0; i < kmer_graph.size(); i++) {
            current_edge = kmer_graph[i];
            current_contig = new ProposedContig(current_edge);
            
            Read* end = current_edge->ext->source;
            Read* start = current_edge->src->source;
            end->usedNTimes++;
            start->usedNTimes++;
            while (true) {
                KMerEdge* best_forward = nullptr;
                int64_t max_weight = 0;
                for (std::vector<KMerEdge*>::size_type i = 0; i < end->kmerEdgesForward.size(); i++) {
                    if (end->kmerEdgesForward[i]->weight > max_weight) {
                        if (end != end->kmerEdgesForward[i]->ext->source) {
                            if (best_forward->ext->source->usedNTimes == 0) {
                                best_forward = end->kmerEdgesForward[i];
                                max_weight = end->kmerEdgesForward[i]->weight;
                            }
                        }
                    }
                }
                if (best_forward != nullptr) {
                    current_contig->addKmerForward(best_forward);
                    end = best_forward->ext->source;
                    end->usedNTimes++;
                }
                KMerEdge* best_backward = nullptr;
                max_weight = 0;
                for (std::vector<KMerEdge*>::size_type i = 0; i < start->kmerEdgesBackward.size(); i++) {
                    if (start->kmerEdgesBackward[i]->weight > max_weight) {
                        if (start != start->kmerEdgesBackward[i]->src->source) {
                            if (best_forward->src->source->usedNTimes == 0) {
                                best_backward = start->kmerEdgesBackward[i];
                                max_weight = start->kmerEdgesBackward[i]->weight;
                            }
                        }
                    }
                }
                if (best_backward != nullptr) {
                    current_contig->addKmerBackward(best_backward);
                    start = best_backward->src->source;
                    start->usedNTimes++;
                }
                //KMer* start = current_contig->kmers.back();
                if (best_forward == nullptr && best_backward == nullptr) {
                    break;
                }
            }
            //printf("EXPORT AT LENGTH %d\n", current_contig->length);
            generated_contig = current_contig->exportContig();
            //printf("> %s\n", generated_contig->sequence);
            
            char name[1024];
            sprintf(name, ">contig-%ld len=%d ofnreads=%ld", i, current_contig->length, current_contig->reads->size());
            output_stream << std::string(name) << "\n";
            output_stream << std::string(generated_contig->sequence) << "\n";
            //fprintf(output_stream, "> contig-%d-len-%d\n", i, current_contig->length);
            //fprintf(output_stream, "%s\n", generated_contig->sequence)
            total_contigs++;
            
            prog++;
            if (prog > (1024 * 16)) {
                printf("\33[2K\r");
                printf("%.3f", (((double) i) / ((double) kmer_graph.size())) * 100.0);
                std::cout << "%" << std::flush;
                prog = 0;
            }
        }
        
        /* for (std::vector<KMerEdge*>::size_type i = 0; i < kmer_graph.size(); i++) {
            
            while (true) {
                //kmer_graph[i]
                break;
            }
            
            prog++;
            if (prog > (1024 * 256)) {
                printf("\33[2K\r");
                printf("%.3f", (((double) i) / ((double) kmer_graph.size())) * 100.0);
                std::cout << "%" << std::flush;
                prog = 0;
            }
        } */
        
        
        
        
        printf("\rMade %ld contigs\n", total_contigs);
        // process above
        
        output_stream << "\n";
        input_stream.close();
    }
    
    if (kontig_mode == 2) {
        
    }
        
    return 0;
}
