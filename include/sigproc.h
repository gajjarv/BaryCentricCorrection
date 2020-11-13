/* sigproc.h: Automatically generated include file for sigproc-4.3 */
#ifdef __cplusplus
extern "C" {
#endif
#include "polyco.h"
#include "epn.h"
#include "version.h"
FILE *open_file(char *filename, char *descriptor) ;
char *backend_name (int machine_id) ;
char *data_category (int data_type) ;
char *headername (char *filename) ;
char *telescope_name (int telescope_id) ;
char tempo_site(int telescope_id) ;
double *bandfactors(int nchans) ;
double *chan_freqs(double fmid, double foff, int nchans, int wapp_off) ;
double deg2dms(double angle) ;
double dmdelay(double f1, double f2, double dm) ;
double ffreq(double tsamp, int npf, int fold, int k) ;
double h2hms(double hours) ;
double inv_cerf(double input) ;
double mjd(int year, int month, int day) ;
double polyco_period(double mjd, struct POLYCO polyco) ;
double polyco_phase(double mjd, struct POLYCO polyco) ;
double pspm_tstart(unsigned long scan_num, char *start_time, double tick_offset, double *mjdobs);
double slaDranrm ( double angle ) ;
double slaGmst ( double ut1 ) ;
double wappcorrect(double mjd) ;
float  *fold_data() ;
float *pspm_prof(FILE *input, int nbins, int nchans, int *table);
float *wapp_prof(int nbins, int nchans, int nifs, int np1, int np2) ;
float flat(float min, float max, long *seed) ;
float gasdev(long *idum) ;
float gauss(long *seed, float mean, float sigma) ;
float nrran0(long *idum) ;
float nrran1(long *idum) ;
float nrran2(long *idum) ;
float nrselect(unsigned long k, unsigned long n, float arr[]) ;
float vmax(float *vec, float n) ;
float vmin(float *vec, float n) ;
int  *bpp_chans(double bw, int mb_start_addr, int mb_end_addr, int mb_start_brd, int mb_end_brd, int *cb_id, double *aib_los, float *dfb_sram_freqs, double rf_lo) ;
int *dmshift(double f1, double df, int nchans, int nbands, double dm, double refrf, double tsamp, double frequency_table[]) ;
int *ignored_channels(char *filename, int nchans) ;
int *pspm_chans(int nchans) ;
int aoftm_read_header(char *filename) ;
int big_endian() ;
int bit(int bitindex, unsigned char byte) ;
int fbin(double tsamp, int npf, int fold, double freq) ;
int file_exists(char *filename) ;
int help_required(char *string) ;
int little_endian() ;
int move_to_keyword(FILE *inputfile, char *keyword) ;
int np2(int n) ;
int process(float time, float start_time, float final_time) ;
int read_block(FILE *input, int nbits, float *block, int nread) ;
int read_blocksubzero(FILE *input, int nbits, float *block, int nread) ;
int read_header(FILE *inputfile) ;
int read_polycoset(FILE *polycofile, struct POLYCO *polyco) ;
int ssm(void) ;
int strings_equal (char *string1, char *string2) ;
int sumchar(unsigned char c) ;
int typeof_inputdata(FILE *fptr, char *filename) ;
int typeof_inputdata(FILE *fptr, char *filename) ;
int vanvleck3lev(double *rho,int npts) ;
long long nsamples(char *filename,int headersize, int nbits, int nifs, int nchans) ;
long long sizeof_file(char name[]) ;
long startseed(void) ;
unsigned char charof2ints (int i, int j) ;
unsigned char charof4ints (int i, int j, int k, int l) ;
void add_channels(float *data, int nsamples, int nadd) ;
void add_samples(float *data, int nifs, int nchans, int nadd) ;
void angle_split(double angle, int *dd, int *mm, double *ss) ;
void aoftm2fb(char *filename, FILE *output) ;
void bandpass_help() ;
void barycentre_help() ;
void blanker_help() ;
void bpp2fb(FILE *input, FILE *output) ;
void cal(double djm, int *year, int *month, int *day) ;
void cel2gal(int rah, int ram, double ras,int decd,int decm,double decs,double *glon,double *glat) ;
void char2fourints (unsigned char c, int *i, int *j, int *k, int *l) ;
void char2ints (unsigned char c, int *i, int *j) ;
void close_log() ;
void create_bro(int nfft) ;
void decimate_data(FILE *input, FILE *output) ;
void decimate_header() ;
void decimate_help() ;
void dedisperse_data(FILE *input, FILE *output) ;
void dedisperse_header() ;
void dedisperse_help() ;
void denorm_prof(float *prof, float *cnt, int nbins, int nifs, int nchans) ;
void depolyco_help() ;
void dice_help() ;
void eight_bit_reorder(unsigned short *spectrum, int nfft) ;
void eraseDP( char *cbuf )	;
void error_message(char *message) ;
void fake_help() ;
void filterbank_header(FILE *outptr) ;
void filterbank_help() ;
void float2char(float *f, int n, float min, float max, unsigned char *c) ;
void float2four(float *f, int n, float min, float max, unsigned char *c) ;
void float2int(float *f, int n, int b, float min, float max, int *i) ;
void float2one(float *f, int n, float min, float max, unsigned char *c) ;
void float2short(float *f, int n, float min, float max, unsigned short *s) ;
void float2two(float *f, int n, float min, float max, unsigned char *c) ;
void fold_header() ;
void fold_help() ;
void four1(float data[], unsigned long nn, int isign) ;
void fshift_prof(float *profile, int nbins, double fshift) ;
void get_nearest_polyco(char *filename, double mjd, struct POLYCO *polyco) ;
void giant_help() ;
void gmrt2fb(FILE *input, FILE *output) ;
void header_help() ;
void headeredit_help() ;
void indexx(unsigned long n, float arr[], unsigned long indx[]) ;
void int2float(int *i, int n, int b, float min, float max, float *f) ;
void machine2prf(FILE *input, FILE *output) ;
void norm_prof(float *prof, float *cnt, int nbins, int nifs, int nchans) ;
void nrsort(unsigned long n, float arr[]) ;
void ooty2fb(FILE *input, FILE *output) ;
void open_log(char *filename) ;
void phcalc(double mjd0,double mjd1,double *phase,double *psrfreq,double *rphase,double *psr_f0,double *poly_tmid,double *coeff,int *num_coeff) ;
void print_version(char *program, char *argument) ;
void prof_adds(float *profile, int *nbins, int nchans, int nifs, int nadd);
void prof_cent(float *profile, int nbins, int nchans) ;
void prof_ddis(float *profile, int nbins, int nchans, int nbands, int nifs, double *chanfreq,  double period, double dm, double reference_frequency, float jyf1, float jyf2) ;
void prof_sbas(char *srcname, float *profile, int nbins, int nchans, int nifs) ;
void prof_sumc(float *profile, int nbins, int nbands, int *nchans, int nifs, int *ignore) ;
void prof_sumifs(float *profile, int nbins, int nchans, int *nifs) ;
void profile_help() ;
void pspm2fb(FILE *input, FILE *output) ;
void pspm_decode(int *rawdata, float *tmparray) ;
void psrfits2fb(FILE *input, FILE *output,char *filename) ;
void pulsar2k2fb(FILE *input, FILE *output) ;
void putd(double d) ;
void putf(float f) ;
void puti(int i) ;
void putl(long l) ;
void putld(long double d) ;
void putu(unsigned long l) ;
void read_aoscan(unsigned long aoscan, int *day, int *year, int *scan) ;
void reader_help() ;
void realft(float data[], unsigned long n, int isign) ;
void rfft(int nfft, double *org, double *fft) ;
void scale_prof(float *profile, int nbins, unsigned long *iprofile, float *scale, float *offset) ;
void scamp2fb(FILE *input, FILE *output) ;
void send_char(char *name, char integer) ;
void send_char_f(char *name, char integer, FILE *output) ;
void send_coords(double raj, double dej, double az, double za) ;
void send_coords_f(double raj, double dej, double az, double za, FILE *output) ;
void send_double (char *name, double double_precision) ;
void send_double_f(char *name, double double_precision, FILE *output) ;
void send_float(char *name,float floating_point) ;
void send_float_f(char *name,float floating_point, FILE *output) ;
void send_int(char *name, int integer) ;
void send_int_f(char *name, int integer, FILE *output) ;
void send_long(char *name, long integer) ;
void send_long_f(char *name, long integer, FILE *output) ;
void send_string(char *string) ;
void send_string_f(char *string, FILE *output) ;
void shift_prof(float *profile, int nbins, int ishift) ;
void slaCaldj ( int iy, int im, int id, double *djm, int *j ) ;
void slaCalyd(int iy, int im, int id, int *ny, int *nd, int *j );
void slaCldj ( int iy, int im, int id, double *djm, int *j ) ;
void slaDjcal ( int ndp, double djm, int iymdf[4], int *j ) ;
void subcal(float *profile, int nbins) ;
void submean(float *profile, int nbins) ;
void submedian(float *profile, int nbins) ;
void sum_fil_help() ;
void swap_double( double *pd ) ;
void swap_float( float *pf ) ;
void swap_int( int *pi ) ;
void swap_long( long *pi ) ;
void swap_longlong( long long *pl ) ;
void swap_short( unsigned short *ps ) ;
void swap_ulong( unsigned long *pi ) ;
void tune_help() ;
void two_bit_reorder(unsigned short *spectrum, int nfft) ;
void update_log(char *string) ;
void uttime(double mjd, int *hh, int *mm, float *ss) ;
void vanvleck9lev(double *rho,int npts) ;
void wapp2fb(FILE *input, FILE *output) ;
void write_dedisp(float *dedisp, int nsout, int nifs, int nbands, float *offset, FILE *output);
void write_epn(FILE *fptr, struct EPN epn) ;
void write_epn_header(FILE *fptr, struct EPN epn) ;
void write_epn_subheader(FILE *fptr, struct EPN epn) ;
void write_profiles(float *prof,int nbins, int nchan, int nifs, FILE *out);
#ifdef __cplusplus
}
#endif
