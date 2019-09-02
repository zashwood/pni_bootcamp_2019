Notes from the hands-on lab running slurm jobs on the spock cluster:

Some prerequisites:
- PNI account so can ssh into scotty interactive server and Spock cluster
- Know how to get around in the Terminal.app on Mac and Putty or other ssh terminal program on PC
  (For PC I recommend installing the official Github command-line client program which installs a bash shell for PC)
- To find Terminal program on Mac, easiest way is to hit Command-Space and type "Terminal" then enter
- Then can add it to your dock by right-clicking on its icon and choosing "Keep in Dock"
- You can play around in Terminal learning to cd, cp, mv, rm, etc -- a great resource is tldp.org bash section
   -- same with the github terminal client on windows, I use it all the time at home. Is a full bash shell.
- xquartz.macosforge.org is where you can get the X11 for Mac environment free. That will allow you to login to
scotty using "ssh -Y" not just "ssh" , so you can see matlab plots or python plots when running those on scotty.
   - we tested that X11 was working via the "xterm" command which brings up a graphical terminal window (mouse-enabled)
   - on PC, there is Xming and something better I couldn't remember. I think it is Mobaxterm from mobaxterm.mobatek.net
   - googling "better than xming" also found VcXsrv which I think might be best sourceforge.net/projects/vcxsrv


-- Interactive use (scotty)

First we went over interactive jobs, which is best done on scotty:

ssh -Y netid@scotty.pni.princeton.edu

You end up in your home directory when you login, which is /usr/people/netid or ~netid as a shortcut (substitute yours).

See https://npcdocs.princeton.edu for more documents on how to access your files from your computer, and lots of other info. And there's always the pnihelp@princeton.edu email for any questions.

--- Screen, resumable interactive sessions.

If you want to be able to close your laptop and walk away, there's the screen command

screen

Anything you do after typing screen will be resumable later; you can see all the screen sessions you have after logging back into scotty via:

screen -list

or just result the only one if you, like me, just have one screen session at most :

screen -r

There's a newer variant called tmux which has more advanced features.

---- Actual job running!

You do this on spock

ssh netid@spock.pni.princeton.edu

Though you can submit jobs from scotty to spock as well.

So there were two example scripts, provided here, ex.sh and ex_arr.sh for one-shot and array jobs, respectively. A one-shot job would be a job that you run once, so it won't run that much faster than it would on your own computer, but you don't have to keep your computer running and can even get email when it's done.

--- the lines we added to get email when the job starts and ends

It isn't in the ex* files because they wouldn't run with "netid" in there , which you will want to replace with yours. Putting in scripts that won't run as-is is annoying I thought. Anyway the lines are:

#SBATCH --mail-type ALL
#SBATCH --mail-user netid@princeton.edu

-- How to actually run the job?

Use the sbatch command:

sbatch ex.sh

-- On the minimal nature of ex.sh and ex_arr.sh

These are the minimum you need to run a job. They only run for 10 minutes tops! After 10 minutes they are killed. It is ok to run for less than 10 (the examples certainly do; probably they complete in only a second or so). So you'd want to at least increase that number in the "SBATCH -t 10" line. They also default to using 1 cpu core and 10G of memory. On npcdocs.princeton.edu, see the Spock* pages for details on other common options. Or google SLURM for an in-depth bit of info, the official docs are at slurm.schedmd.com. We are running version 18.08 as of August 2019 at the PNI and the SLURM docs are about version 19 as of that time. Most info will be the same.

--- Array jobs: ex_arr.sh

The ex_arr.sh file adds another SBATCH line:
#SBATCH --array 1-2

This just says to actually run ex_arr.sh twice, once with the taskid set to 1 , another with it set to 2 (it could run 2 before 1 so cannot depend on the order in your code). You'll need to check the value of SLURM_ARRAY_TASK_ID environment variable in your script, or pass that value into your script, to do something contingent on its value. Then you can really get stuff done because you could have hundreds or thousands of concurrent jobs running! Each can check its number, load a particular data set, or explore its own part of parameter space. Now you are clusterin'

Anyway sorry for the quick and informal nature of this intro but again you can get more info at npcdocs.princeton.edu, the princeton.edu website (including info on princeton's matlab licenses and how to install, as someone asked) and there's always pnihelp@princeton.edu. And I'm in B14B down by the glass double door near the parking lot. And there are cluster help sessions for the entire campus on Tuesdays 10:30-11:30 in Lewis Library 245 and 2-3 on Thursdays same place. I volunteer on the tuesdays so might see you there, or by appt in my office. Have fun!

-Ben Singer, research applications specialist , PNI  <bdsinger@princeton.edu>
