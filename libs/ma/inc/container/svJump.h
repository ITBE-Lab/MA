#pragma once

class SvJump
{
  public:
    enum SeedOrientation
    {
        forwardToForward,
        reverseToReverse,
        forwardToReverse,
        reverseToForward
    }; // enum
    nucSeqIndex uiFrom;
    nucSeqIndex uiTo;
    nucSeqIndex uiQueryDistance;
     xSeedOrientation;
}; // class