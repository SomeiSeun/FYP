# FYP
Code and logs from my thesis can be found in this repo.

Hello! If you are not the next person to do this FYP, uhhh imagine me talking to you makes sense, or just see ths stuff. Anyways, This will be a really interesting and fun FYP and I hope I've left you with an interesting starting point for your FYP.

So you will find here these folders:

- FYP-notes -> Essentially my logbook during the FYP, not particularly important but you can import it into Obsidian as a repository and see what I noted down. I mainly recommend using obsidian as a logbook type tool, it really helps with keeping track of each aspect.
- GMAT-FILES -> Relevant GMAT files and matlab SPAD generation files for your use.
- MyWork -> All of my computational work from the thesis is in this folder.


To see what each code does, you can find it in the appendix of my report. I will highlight some files to look at and use:
- TorqueOnSailOffset: With this function you can compute the torque on any size of sail and it will calculate the center of pressure based on orientation and position of the sail. You can adapt this code to account for other axis rotations as well, potentially reworking the rotation implementation using RotZYZ function.
- PowerGenCon: Code from which I performed the power subsystem sizing. It wasn't an optimisation, more trial and error and intiution.
- HubSizing: Code to calculate the clearance distance between a stowed sail and the inner wall. You can adjust parameters to the appropriate size you are testing on.

I am sure Dr Knoll will debrief you on important findings and considerations going forward. He had a really interesting idea about a custom device to adjust the CoG placement to correct any induced rotations. You might have to rework the equations of the resultant acceleration that I calculated in EulerEq12 (sorry about that..)

In the FutureStudyConsideration folder inside the Code folder, you will find a lil file detailing some additional thoughts I had for going forward in the project that I didn't fully expand upon in my thesis itself.
