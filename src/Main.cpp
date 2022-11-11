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

int main(int argc, char** argv) {
    std::cout << "kontig FASTQ contig generator\n\n";
    char* input_file = NULL;
    char* output_file = NULL;

    bool stdout_switch = true;
    
    for (int i = 1; i < argc; ++i) {
        char* arg = argv[i];
        if (strcmp(arg,"-i") == 0) {
            i++;
            input_file = argv[i];
            continue;
        }
        if (strcmp(arg,"-o") == 0) {
            i++;
            output_file = argv[i];
            stdout_switch = false;
            continue;
        }
        if (strcmp(arg,"-stdout") == 0) {
            stdout_switch = true;
            continue;
        }
        //bootVmImmediate = arg;
    }

    if (input_file == NULL) {
        std::cout << "Need to specify input file with '-i'.\n";
        exit(1);
    }
    
    if (stdout_switch) {
        std::cout << input_file << " -> stdout\n";
    } else {
        std::cout << input_file << " -> " << output_file << "\n";
    }

    std::ifstream input_stream(input_file);
    
    std::unique_ptr<std::ostream> output_stream = std::make_unique<std::ostream>(std::cout.rdbuf());
    
    if (!stdout_switch) {
        output_stream = std::make_unique<std::ofstream>(output_file);
    }
    
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
        
        input_stream.read(input_buffer, 1024);
        chars_read = input_stream.gcount();
        still_readable = chars_read > 0;
    } while (still_readable);
    
    //(*output_stream).write(total_output, total_output_length);
    
    free(input_buffer);
    //exit(0);
    
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
    printf("Read %d characters from FASTQ file\n", input_position);
    
    // process below
    
    KMerSet** kmer_set_set = (KMerSet**) malloc(sizeof(KMerSet*) * reads.size());
    int total_kmers = 0;
    
    int kmer_length = 128;
    
    for (int i = 0; i < reads.size(); i++) {
        //printf("%s : %d", read->name, read->length);
        struct KMerSet* kms = reads[i]->generateKmers(kmer_length); // 4^n  = 2^2n
        total_kmers += kms->length;
        kmer_set_set[i] = kms;
        //printf("%s : %d", read->sequence, kms->length);
        //*output_stream << "\n";
    }
    
    printf("Generated %d kmers\n", total_kmers);
    
    KMer** kmers = (KMer**) malloc(sizeof(KMer*) * total_kmers);
    
    int current_kmer = 0;
    
    for (int i = 0; i < reads.size(); i++) {
        for (int j = 0; j < kmer_set_set[i]->length; j++) {
            
            kmers[current_kmer] = kmer_set_set[i]->kmers[j];
            current_kmer++;
        }
    }
    
    free(kmer_set_set);
    
    //for (int i = 0; i < total_kmers; i++) {
        //printf("%.16s\n", kmers[i]->sequence);
    //}
    
    printf("%d K-Mers about to be sorted\n", total_kmers);
    
    uint64_t prog = 0;
    
    // KMers prepared
    // Now sort KMers
    
    KMerSortingTreeNode* root_kmer_node = new KMerSortingTreeNode();
    root_kmer_node->high = nullptr;
    root_kmer_node->low = nullptr;
    root_kmer_node->values = nullptr;
    
    
    
    for (int i = 0; i < total_kmers; i++) {
        KMerSortingTreeNode* cnode = root_kmer_node;
        for (int j = 0; j < 64; j++) {
            if ((kmers[i]->quickRef >> j) & 0b1) {
                if (cnode->high == nullptr) {
                    cnode->high = new KMerSortingTreeNode();
                    cnode->high->high = nullptr;
                    cnode->high->low = nullptr;
                    cnode->high->values = nullptr;
                }
                cnode = cnode->high;
            } else {
                if (cnode->low == nullptr) {
                    cnode->low = new KMerSortingTreeNode();
                    cnode->low->high = nullptr;
                    cnode->low->low = nullptr;
                    cnode->low->values = nullptr;
                }
                cnode = cnode->low;
            }
        }
        KMerLinkedList* ll = cnode->values;
        if (ll == nullptr) {
            ll = new KMerLinkedList();
            ll->next = nullptr;
            ll->value = kmers[i];
            cnode->values = ll;
        } else {
            while (ll->next != nullptr) {
                ll = ll->next;
            }
            ll->next = new KMerLinkedList();
            ll->next->next = nullptr;
            ll->next->value = kmers[i];
        }
        
        prog++;
        if (prog > (1024 * 64)) {
            printf("\33[2K\r");
            printf("%.3f", (((double) i) / ((double) total_kmers)) * 100.0);
            std::cout << "%" << std::flush;
            //fflush();
            prog = 0;
        }
    }
    
    printf("\rK-Mers loaded into tree\n");
    
    // Now connect KMers

    std::vector<KMerEdge*> kmer_graph;
    
    prog = 0;
    
    printf("Connecting K-Mers\n");
    
    uint64_t total_connections = 0;
    
    for (int i = 0; i < reads.size(); i++) {
        KMer* head_kmer = reads[i]->kmers->kmers[0];
        //KMer* tail_kmer = reads[i]->kmers->kmers[reads[i]->kmers->length];
        //head_kmer->quickRef;
        KMerSortingTreeNode* cnode = root_kmer_node;
        for (int j = 0; j < 64; j++) {
            if ((head_kmer->quickRef >> j) & 0b1) {
                cnode = cnode->high;
            } else {
                cnode = cnode->low;
            }
        }
        //if (cnode->values == nullptr) {
        //    printf("uwu");
        //}
        KMerLinkedList* ll = cnode->values;
        while (ll != nullptr) {
            KMer* candidate_kmer = ll->value;
            if (head_kmer->source != candidate_kmer->source) {
                //printf("0x%08lx == 0x%08lx\n", kmers[i]->quickRef, kmers[j]->quickRef);
                if (head_kmer->isEqual(candidate_kmer)) {
                    KMerEdge* kmer_node = new KMerEdge();
                    kmer_node->src = candidate_kmer;
                    kmer_node->ext = head_kmer;
                    //kmer_node->weight = head_kmer->getMatchQuality(kmers[j]);
                    kmer_graph.push_back(kmer_node);
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
            ll = ll->next;
        }
    }
    
    printf("\rMade %ld connections\n", total_connections);
    
    // Connections made, time to make proposed contigs via greedy algorithm
    
    
    
    // process above
    
    *output_stream << "\n";
    
    
    
    input_stream.close();
    
    return 0;
}
