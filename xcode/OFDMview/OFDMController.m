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
	[model adjust :([freqadjust doubleValue] / 10.0)];

	[nsamples setStringValue:@"Running FFT..."];
	[self rerunFFT :sender];
	[nsamples setIntValue:[model getSampleCount]];
}

- (IBAction) rerunFFT : (id) sender {
	[lowerspec setImage: [[NSImage alloc] initWithSize: NSMakeSize(150, 150)]];
	[upperspec setImage: [[NSImage alloc] initWithSize: NSMakeSize(150, 150)]];

	[model runFFT :[offset intValue] :[lowerspec image] :[upperspec image]];
	[self rerenderConstellation:sender];
}

- (IBAction) rerenderConstellation : (id) sender {
	[iview setImage: [[NSImage alloc] initWithSize: NSMakeSize(150, 150)]];
	[model constellationIter :self :[carrier intValue] :[iview image]];
}

@end
