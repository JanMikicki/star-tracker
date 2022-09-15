# star-tracker
Emulating star-tracker grid identification algorithm.

Star trackers try to determine the orientation (or attitude) of the spacecraft with respect to the stars.

The algorithm receives a random rotated tile from the sky map (lower left image).
Orange line represents star-tracker's local top/north.

The algorithm guesses the stars ids (lower right) in this view by comparing it against the entire sky map.

![image](https://user-images.githubusercontent.com/13591149/190478483-f4106a18-e0d4-4dd2-a9ad-882dc7292604.png)
