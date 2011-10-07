// OFDMModel definition
// OFDMview

#import <Cocoa/Cocoa.h>
#include <fftw3.h>

@interface OFDMModel : NSObject {
	double *samples;
	int nsamples;
	fftw_complex *symbols;
	int nsymbols;
}

- (int)loadData :(char*)filename;
- (void)adjust : (double)frq;
- (int)getSampleCount;
- (void)runFFT : (int)startsamp :(id)lower :(id)upper;
- (void)constellationIter :(id)sender :(int)carrier :(id)constellation;

@end
