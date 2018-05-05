import numpy as np

from sklearn import linear_model, datasets
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.layouts import row, column, layout
from bokeh.palettes import d3
from bokeh.models import LinearAxis, Range1d, LogColorMapper, FixedTicker, BasicTicker, Grid
from bokeh.models.formatters import FuncTickFormatter
from bokeh.layouts import gridplot
from bokeh.io import save


print(list(range(0, 10, 1)))

n_samples = 1000
n_outliers = 50




X, y, coef = datasets.make_regression(n_samples=n_samples, n_features=1,
                                      n_informative=1, noise=10,
                                      coef=True, random_state=0)

# Add outlier data
np.random.seed(0)
X[:n_outliers] = 3 + 0.5 * np.random.normal(size=(n_outliers, 1))
y[:n_outliers] = -3 + 10 * np.random.normal(size=n_outliers)

X[n_outliers:n_outliers*2] = 0 + 0.5 * np.random.normal(size=(n_outliers, 1))
y[n_outliers:n_outliers*2] = 200 + 10 * np.random.normal(size=n_outliers)



test_x = [0, 1, 2, 3]
test_y = [0, 1, 2, 3]

X = np.array(test_x)[:, np.newaxis]
y = np.array(test_y)


# Fit line using all data
lr = linear_model.LinearRegression()
lr.fit(X, y)

# Robustly fit linear model with RANSAC algorithm
ransac = linear_model.RANSACRegressor()
ransac.fit(X, y)
inlier_mask = ransac.inlier_mask_
outlier_mask = np.logical_not(inlier_mask)

# Predict data of estimated models
line_x = np.arange(X.min(), X.max())
line_X = line_x[:, np.newaxis]
line_y = lr.predict(line_X)
line_y_ransac = ransac.predict(line_X)

# Compare estimated coefficients
print("Estimated coefficients (true, linear regression, RANSAC):")
print(coef, lr.coef_, ransac.estimator_.coef_)

lw = 2
#line_Xr = range(X.min(), X.max())[:, np.newaxis]
print("===========")
plot = figure(width=800)
print(np.newaxis)
print(X.flatten())
plot.x(X.flatten(), y, color="blue")
plot.line(line_x, line_y, color="green", line_width=4)
plot.line(line_x, line_y_ransac, color="red", line_width=4)

show(plot)


#plt.scatter(X[inlier_mask], y[inlier_mask], color='yellowgreen', marker='.',
#            label='Inliers')
#plt.scatter(X[outlier_mask], y[outlier_mask], color='gold', marker='.',
#            label='Outliers')
#plt.plot(line_X, line_y, color='navy', linewidth=lw, label='Linear regressor')
#plt.plot(line_X, line_y_ransac, color='cornflowerblue', linewidth=lw,
#         label='RANSAC regressor')
#plt.legend(loc='lower right')
#plt.xlabel("Input")
#plt.ylabel("Response")
#plt.show()