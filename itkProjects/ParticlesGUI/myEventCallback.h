//
//  myEventCallback.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/17/12.
//
//

#ifndef ParticlesGUI_myEventCallback_h
#define ParticlesGUI_myEventCallback_h

class EventCallback {
public:
    virtual void EventRaised(int eventId, int eventCode, const void* src = 0, void* data = 0) = 0;
};

#endif
