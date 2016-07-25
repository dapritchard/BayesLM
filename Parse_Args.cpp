
#include <cstdlib>    // strncmp, strlen
#include <cstring>    // strtol, strtod
#include <cerrno>
#include <stdexcept>  // runtime_exception
#include <getopt.h>

#include <iostream>

#include "Parse_Args.h"  // FILENAME_MAXLEN definition

#define MAX_ERR_LEN 256  // Amount of memory allocated for error messages

// See Bayes_LM_Arma.cpp for description of these variables
extern int n;
extern int p;
extern double prop_nonzero;
extern double true_beta_sd;
extern double true_sigma;
extern double sd_x;
extern int nsamp;
extern bool print_stats;
extern bool write_ctime;
extern char ctime_file_loc[];
extern bool write_samples;
extern char samples_file_loc[];
extern char decomp_method;
extern unsigned int seed;

void set_param(int param_idx, char* arg_ptr);




// Argument parsing functions --------------------------------------------------


/* Parse user parameter specifications and set the appropriate variable values.
 * Input of arguments is allowed following the usual UNIX conventions.
 */

void parse_args(int argc, char* argv[]) {

    int c;                      // Argument type info
    int option_idx;             // Option index for long_options array
    char err_msg[MAX_ERR_LEN];  // Memory for error message

    struct option long_options[] = {
	{ "ctime_file_loc",   required_argument, 0, 0 },  // option_idx  0
	{ "decomp_method",    required_argument, 0, 0 },  // option_idx  1
	{ "n",                required_argument, 0, 0 },  // option_idx  2
	{ "nsamp",            required_argument, 0, 0 },  // option_idx  3
	{ "p",                required_argument, 0, 0 },  // option_idx  4
	{ "print_stats",      required_argument, 0, 0 },  // option_idx  5
	{ "prop_nonzero",     required_argument, 0, 0 },  // option_idx  6
	{ "samples_file_loc", required_argument, 0, 0 },  // option_idx  7
	{ "sd_x",             required_argument, 0, 0 },  // option_idx  8
	{ "seed",             required_argument, 0, 0 },  // option_idx  9
	{ "true_beta_sd",     required_argument, 0, 0 },  // option_idx 10
	{ "true_sigma",       required_argument, 0, 0 },  // option_idx 11
	{ "write_ctime",      required_argument, 0, 0 },  // option_idx 12
	{ "write_samples",    required_argument, 0, 0 },  // option_idx 13
	{ 0,                  0,                 0, 0 },  // required last val
    };

    while (1) {

	c = getopt_long_only(argc, argv, ":", long_options, &option_idx);
	if (c == -1) {
	    break;
	}

	switch (c) {
	case 0:
	    set_param(option_idx, optarg);
	    break;
	case '?':
	    strcat(err_msg, "Invalid option ");
	    // 15 is the length of "Invalid option "
	    strncat(err_msg + 15, argv[--optind], MAX_ERR_LEN - 16);
	    // Ensure string is null-terminated
	    err_msg[MAX_ERR_LEN - 1] = '\0';
	    std::cout << "err message is:  " << err_msg << "\n";
	    throw std::runtime_error(err_msg);
	case ':':
	    strcat(err_msg, "Missing argument for ");
	    // 22 is the length of "Missing argument for "
	    strncat(err_msg + 21, argv[--optind], MAX_ERR_LEN - 22);
	    // Ensure string is null-terminated
	    err_msg[MAX_ERR_LEN - 1] = '\0';
	    throw std::runtime_error(err_msg);
	}
    }
}




/* Process argument pointed to by arg_ptr, where param_idx tells us which
 * variable the argument is providing a value for, and set the appropriate
 * global variable to the provided value.
 */

void set_param(int param_idx, char* arg_ptr) {

    // Point to next char after a int / double read (should point to '\0')
    char* endptr;
    // Temporary integer for storage of a read by strtol
    int temp_int;
    // To distinguish success / failure after a call to strtol or strtod
    errno = 0;

    /* Map of param_idx given below:
     *
     *      0: ctime_file_loc
     *      1: decomp_method
     *      2: n
     *      3: nsamp
     *      4: p
     *      5: print_stats
     *      6: prop_nonzero
     *      7: samples_file_loc
     *      8: sd_x
     *      9: seed
     *     10: true_beta_sd
     *     11: true_sigma
     *     12: write_ctime
     *     13: write_samples
     */

    switch (param_idx) {

    case 0:  // ctime_file_loc
	if (strlen(arg_ptr) > FILENAME_MAXLEN - 1) {
	    throw std::runtime_error("ctime_file_loc too long\n");
	}
	strcpy(ctime_file_loc, arg_ptr);
	break;
    case 1:  // decomp_method
	if (strcmp(arg_ptr, "chol") == 0) decomp_method = 'c';
	else if (strcmp(arg_ptr, "eigen") == 0) decomp_method = 'e';
	else throw std::runtime_error("method must be either \"chol\" or \"eigen\"\n");
	break;
    case 2:  // n
	n = strtol(arg_ptr, &endptr, 10);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for n\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow for n\n");
	else if (n <= 0) throw std::runtime_error("n must be > 0\n");
	break;
    case 3:  // nsamp
	nsamp = strtol(arg_ptr, &endptr, 10);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for nsamp\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow for nsamp\n");
	else if (nsamp <= 0) throw std::runtime_error("nsamp must be > 0\n");
	break;
    case 4:  // p
	p = strtol(arg_ptr, &endptr, 10);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for p\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow p\n");
	else if (p <= 0) throw std::runtime_error("p must be > 0\n");
	break;
    case 5:  // print_stats
	if (strcmp(arg_ptr, "true") == 0) print_stats = true;
	else if (strcmp(arg_ptr, "false") == 0) print_stats = false;
	else throw std::runtime_error("print_stats must be either \"true\" or \"false\"\n");
	break;
    case 6:  // prop_nonzero
	prop_nonzero = strtod(arg_ptr, &endptr);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for prop_nonzero\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow for prop_nonzero\n");
	else if (prop_nonzero <= 0) throw std::runtime_error("prop_nonzero must be > 0\n");
	else if (prop_nonzero > 1) throw std::runtime_error("prop_nonzero must be <= 1\n");
	break;
    case 7:  // samples_file_loc
	if (strlen(arg_ptr) > FILENAME_MAXLEN - 1) {
	    throw std::runtime_error("samples_file_loc too long\n");
	}
	strcpy(samples_file_loc, arg_ptr);
	break;
    case 8:  // sd_x
	sd_x = strtod(arg_ptr, &endptr);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for sd_x\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow for sd_x\n");
	else if (sd_x <= 0) throw std::runtime_error("sd_x must be > 0\n");
	break;
    case 9:  // seed
	temp_int = strtol(arg_ptr, &endptr, 10);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for seed\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow seed\n");
	else if (temp_int < 0) throw std::runtime_error("seed must be >= 0\n");
	seed = temp_int;
	break;
    case 10:  // true_beta_sd
	true_beta_sd = strtod(arg_ptr, &endptr);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for true_beta_sd\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow for true_beta_sd\n");
	else if (true_beta_sd <= 0) throw std::runtime_error("true_beta_sd must be > 0\n");
	break;
    case 11:  // true_sigma
	true_sigma = strtod(arg_ptr, &endptr);
	if (*endptr != '\0') throw std::runtime_error("Invalid argument for true_sigma\n");
	else if (errno != 0) throw std::runtime_error("Underflow / overflow for true_sigma\n");
	else if (true_sigma <= 0) throw std::runtime_error("true_sigma must be > 0\n");
	break;
    case 12:  // write_ctime
	if (strcmp(arg_ptr, "true") == 0) write_ctime = true;
	else if (strcmp(arg_ptr, "false") == 0) write_ctime = false;
	else throw std::runtime_error("write_ctime must be either \"true\" or \"false\"\n");
	break;
    case 13:  // write_samples
	if (strcmp(arg_ptr, "true") == 0) write_samples = true;
	else if (strcmp(arg_ptr, "false") == 0) write_samples = false;
	else throw std::runtime_error("write_samples must be either \"true\" or \"false\"\n");
	break;

    } // end parse argument switch

}




