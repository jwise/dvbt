//
//  OFDMController.m
//  OFDMview
//
//  Created by Joshua Wise on 7/26/09.
//  Copyright 2009 __MyCompanyName__. All rights reserved.
//

#import "OFDMController.h"


@implementation OFDMController

- (id) init {
	self = [super init];
	model = [[OFDMModel alloc] init];
	
	return self;
}

- (IBAction) loadFile : (id) sender {
	

	[nsamples setStringValue:@"Loading data..."];
	[model loadData:"/Users/joshua/projects/dvbt/dvbt.mixed.raw"];
	
	[nsamples setStringValue :@"Frequency adjusting..."];
	[model adjust :[freqadjust intValue]];

	[nsamples setStringValue:@"Running FFT..."];
	[model runFFT :[offset intValue]];
	[nsamples setIntValue:[model getSampleCount]];
	[self rerenderConstellation:sender];
}

- (IBAction) rerunFFT : (id) sender {
	[model runFFT :[offset intValue]];
	[self rerenderConstellation:sender];
}

- (IBAction) rerenderConstellation : (id) sender {
	[iview setImage: [[NSImage alloc] initWithSize: NSMakeSize(150, 150)]];
	
	[[iview image] lockFocus];
	[model constellationIter :self :[carrier intValue]];
	[[iview image] unlockFocus];
}

@end
