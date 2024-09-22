This program takes a list of numbers and uses its moments (mean, variance, skew, etc.)
    to guess which probability distribution the data follows.
It assumes the sample moments are approximately the population moments,
    so don't use it on small data sets.
I'm not sure how the results change in response to small changes in the data.
The input file should be just a list of floats separated by newlines.

The data is given a moment score for each distribution:
    sum((sample moment - moment assuming distribution)^2)/(# of independent moments)
    A moment score closer to 0 indicates the relationship between the moments
    lines up more closely with that distribution.
