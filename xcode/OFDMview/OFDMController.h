//
//  OFDMController.h
//  OFDMview
//
//  Created by Joshua Wise on 7/26/09.
//  Copyright 2009 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "OFDMModel.h"

@interface OFDMController : NSObject {
	IBOutlet NSTextField *nsamples;
	IBOutlet NSImageView *iview;
	IBOutlet NSImageView *lowerspec;
	IBOutlet NSImageView *upperspec;
	IBOutlet NSStepper *carrier;
	IBOutlet NSSlider *offset;
	IBOutlet NSSlider *freqadjust;
	
	OFDMModel *model;
}

- (IBAction)loadFile : (id)sender;
- (IBAction)rerunFFT : (id)sender;
- (IBAction)rerenderConstellation : (id)sender;

@end
