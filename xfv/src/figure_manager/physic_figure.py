# -*- coding: utf-8 -*-
"""
Class to define a figure
"""
from os import sep
import matplotlib.pyplot as plt


class PhysicFigure:
    """
    Figure
    """
    def __init__(self, X, Y, xlabel="X", ylabel="Y", titre="titre", interface_id=0, save_path=None):
        self._fig, self._ax = plt.subplots()
        self._line_target, = self._ax.plot(X, Y, '-+', color='blue')
        self._line_projectile, = self._ax.plot(X, Y, '-+', color='purple')
        self._ax.set_xlabel(xlabel)
        self._ax.set_ylabel(ylabel)
        self._ax.set_title(titre)
        self._fig_number = 1
        self._title = titre
        self._fig.canvas.draw()
        self._save_path = save_path
        self.interface_id = interface_id  # to separate target and projectile
        plt.show(block=False)

    def set_y_limit(self, val_min=0., val_max=1.0):
        """
        Fixation des limites en y
        """
        self._ax.set_ylim([val_min, val_max])

    def set_x_limit(self, val_min=0., val_max=1.0):
        """
        Fixation des limites en x
        """
        self._ax.set_xlim([val_min, val_max])

    def update(self, abscissa=None, ordinate=None, title_comp=None):
        """
        Update image to get animation
        And save the image if a path is given
        """
        if abscissa is not None:
            self._line_target.set_xdata(abscissa[self.interface_id:])
            self._line_projectile.set_xdata(abscissa[:self.interface_id])
        if ordinate is not None:
            self._line_target.set_ydata(ordinate[self.interface_id:])
            self._line_projectile.set_ydata(ordinate[:self.interface_id])
        if title_comp is not None:
            self._ax.set_title(self._title + ' ' + title_comp)

        self._ax.draw_artist(self._ax.patch)
        self._ax.draw_artist(self._line_target)
        self._ax.draw_artist(self._line_projectile)
        self._fig.canvas.update()
        self._fig.canvas.flush_events()

        if self._save_path is not None:
            rac_path = self._save_path + sep + self._title
            fig_path = rac_path + "_{:04d}.png".format(self._fig_number)
            fig_path = fig_path.replace(" ", "_")
            data_path = rac_path + "_{:04d}.dat".format(self._fig_number)
            data_path = data_path.replace(" ", "_")
            self._fig.savefig(fig_path)
            self._fig_number += 1
            with open(data_path, "w") as file_object:
                for x_data, y_data in zip(abscissa, ordinate):
                    file_object.write("{:20.18g}{:s}{:20.18g}\n".format(float(x_data), " ",
                                                                        float(y_data)))
