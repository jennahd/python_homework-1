"""
We all have the problem of too many songs and too little time...
Parse your song lists and you can immediately know how long a song 
you're interested in is and enjoy a youtube video of that song!

All you need to do is provide a list of songs in your home directory!
Each entry should include:
1. song title
2. song artist
3. song duration
"""

import os
import sh
import warnings
import webbrowser

class Song(object):
    """Describes the attributes (re. artist and length) of
    different songs"""
    def __init__(self, title, artist, duration):
        self.title = title
        self.artist = artist
        #Check that duration is okay#
        try:
            self.duration = int(duration)
        except ValueError:
            warnings.warn("The song '%s' is not a song duration!! It will be set to 0." % self.title)
            self.duration = 0
        # Check duration is not negative #
        if self.duration < 0:
            raise Exception("You can't have a negative duration on song...change this please! Cannot play: " + self.title)
    
    def pretty_duration(self):
        seconds = self.duration
        minutes = seconds / 60
        hours = minutes / 60
        return "%02i hours %02i minutes %02i seconds" % (hours, minutes % 60, seconds % 60)
    def play(self):
        url = "https://www.youtube.com/results?search_query="
        final_url = url + self.title.replace(' ','+')
        webbrowser.open(final_url)

location = os.environ['HOME']
in_path = location + '/lulu_mix_16.csv'
if not os.path.exists(in_path):
    raise Exception("No file at '%s'." % in_path)

song_info = [Song(*l.strip('\n').split(',')) for l in open(in_path, 'rb')]


#Test Code
for s in song_info: print s.artist
for s in song_info: print s.pretty_duration()
print sum(s.duration for s in song_info), "seconds in total"
song_info[6].play()
