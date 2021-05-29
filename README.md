# Smallest circumscribing *k*-gon

Implementation of the smallest circumscribing *k*-gon algorithm presented in [ACY85].

Given a convex *n*-gon and a number *k*, this algorithm can find the smallest *k*-gon that circumscribes the *n*-gon in O(*n*<sup>2</sup> log*n* *k*) time.

I mostly use it for [Particle Trimming](http://www.humus.name/index.php?ID=266), but you may find it useful in other places. Usage can be found in <example/main.cpp>.

---

[ACY85] Aggarwal, Alok, Jyun-Sheng Chang, and Chee K. Yap. "Minimum area circumscribing polygons." The Visual Computer 1.2 (1985): 112-117.