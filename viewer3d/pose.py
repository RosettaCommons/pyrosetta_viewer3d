import attr
import collections
import logging

from ipywidgets.widgets import Output
from pyrosetta import Pose
from typing import Iterable, List, NoReturn, Optional, Union

from viewer3d.validators import requires_show


_logger: logging.Logger = logging.getLogger("viewer3d.pose")
out = Output()


@attr.s(kw_only=False, slots=False)
class PoseBase:
    @requires_show
    def add_pose(
        self, pose: Pose, index: Optional[int] = None, update_viewer: bool = True
    ) -> Optional[NoReturn]:
        """
        Add a pose object to the Viewer instance, optionally updating the visualization.

        Args:
            pose: a required `Pose` object representing the pose to add to the viewer.
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None`
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if the pose is not of type `Pose`.
            `AssertionError` if the index is not of type `int`.

        Returns:
            `None`
        """
        assert isinstance(
            pose, Pose
        ), f"Object must be of type `Pose`. Received: {type(pose)}"
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        self.poses[index].append(pose)
        self.pdbstrings[index].append(None)
        if update_viewer:
            model = self.get_last_model(index)
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=False
            )

    @requires_show
    def add_pdbstring(
        self, pdbstring: str, index: Optional[int] = None, update_viewer: bool = True
    ) -> Optional[NoReturn]:
        """
        Add a PDB string object to the Viewer instance, optionally updating the visualization.

        Args:
            pdbstring: a `str` object representing the PDB string to add to the viewer.
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None`
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if PDB string is not of type `str`.
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        assert isinstance(
            pdbstring, str
        ), f"Object must be of type `str`. Received: {type(pdbstring)}"
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        self.poses[index].append(None)
        self.pdbstrings[index].append(pdbstring)
        if update_viewer:
            model = self.get_last_model(index)
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=False
            )

    def remove_pose(
        self,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        """
        Remove a pose object from the Viewer instance, optionally updating the visualization.

        Args:
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None`
            model: an optional `int` object representing the model in the poses or pdbstrings
                index to update. If `None`, then update all models in the poses or pdbstrings
                index to update.
                Default: `None`
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        self.remove_pdbstring(index=index, model=model, update_viewer=update_viewer)

    @requires_show
    def remove_pdbstring(
        self,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        """
        Remove a PDB string object from the Viewer instance, optionally updating the visualization.

        Args:
            index: an optional `int` object representing the pdbstrings index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None`
            model: an optional `int` object representing the model in the poses or pdbstrings
                index to update. If `None`, then update all models in the poses or pdbstrings
                index to update.
                Default: `None`
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        if model is None or model not in set(range(len(self.pdbstrings[index]))):
            model = self.get_last_model(index)
        if len(self.pdbstrings[index]) > 0:
            self.poses[index].pop(model)
            self.pdbstrings[index].pop(model)
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes are empty at index `{index}`."
            )
        if update_viewer:
            self.update_viewer(
                index=index, model=model, add_objects=False, remove_objects=True
            )

    @requires_show
    def update_pose(
        self,
        pose: Pose,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        """
        Update a pose object in the Viewer instance, optionally updating the visualization.

        Args:
            pose: a `Pose` object representing the pose to update in the viewer.
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None`
            model: an optional `int` object representing the model in the poses or pdbstrings
                index to update. If `None`, then update all models in the poses or pdbstrings
                index to update.
                Default: `None`
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if pose is not of type `Pose`.
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        assert isinstance(
            pose, Pose
        ), f"Object must be of type `Pose`. Received: {type(pose)}"
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        if model is None or model not in set(range(len(self.poses[index]))):
            model = 0
        if index in self.poses.keys():
            self.poses[index][model] = pose
            self.pdbstrings[index][model] = None
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )
        if update_viewer:
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=True
            )

    @requires_show
    def update_pdbstring(
        self,
        pdbstring: str,
        index: Optional[int] = None,
        model: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        """
        Update a PDB string object in the Viewer instance, optionally updating the visualization.

        Args:
            pdbstring: a `str` object representing the PDB string to update in the viewer.
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None`
            model: an optional `int` object representing the model in the poses or pdbstrings
                index to update. If `None`, then update all models in the poses or pdbstrings
                index to update.
                Default: `None`
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if PDB string is not of type `str`.
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        if model is None or model not in set(range(len(self.pdbstrings[index]))):
            model = 0
        if index in self.pdbstrings.keys():
            self.poses[index][model] = None
            self.pdbstrings[index][model] = pdbstring
        else:
            raise IndexError(
                f"The 'poses' and 'pdbstrings' attributes do not have index `{index}`."
            )
        if update_viewer:
            self.update_viewer(
                index=index, model=model, add_objects=True, remove_objects=True
            )

    @requires_show
    def update_poses(
        self, poses: List[Pose], index: Optional[int] = None, update_viewer: bool = True
    ) -> Optional[NoReturn]:
        """
        Update pose objects at a certain index in the Viewer instance, optionally
        updating the visualization.

        Args:
            poses: a `list` object of `Pose` objects representing the poses to update in
                the viewer.
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None``
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if poses object is not of type `list`.
            `AssertionError` if pose object is not of type `Pose`.
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        assert isinstance(poses, list), "The 'poses' argument must be of type `list`."
        for pose in poses:
            assert isinstance(
                pose, Pose
            ), f"Object must be of type `Pose`. Received: {type(pose)}"
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        self.poses[index] = poses
        self.pdbstrings[index] = [None] * len(poses)
        if update_viewer:
            self.update_viewer(
                index=index, model=None, add_objects=True, remove_objects=True
            )

    @requires_show
    def update_pdbstrings(
        self,
        pdbstrings: List[str],
        index: Optional[int] = None,
        update_viewer: bool = True,
    ) -> Optional[NoReturn]:
        """
        Update PDB string objects at a certain index in the Viewer instance, optionally
        updating the visualization.

        Args:
            pdbstrings: a `list` object of `str` objects representing the PDB strings to update
                in the viewer.
            index: an optional `int` object representing the poses index to update.
                If `None`, then update all poses in the displayed index.
                Default: `None``
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Raises:
            `AssertionError` if PDB strings object is not of type `list`.
            `AssertionError` if PDB string object is not of type `str`.
            `AssertionError` if index is not of type `int`.

        Returns:
            `None`
        """
        assert isinstance(
            pdbstrings, list
        ), "The 'pdbstrings' argument must be of type `list`."
        for pdbstring in pdbstrings:
            assert isinstance(
                pdbstring, str
            ), f"Object must be of type `str`. Received: {type(pdbstring)}"
        if index is None:
            index = self.get_decoy_widget_index()
        else:
            assert isinstance(index, int), "Index must be of type `int`."
        self.poses[index] = [None] * len(pdbstrings)
        self.pdbstrings[index] = pdbstrings
        if update_viewer:
            self.update_viewer(
                index=index, model=None, add_objects=True, remove_objects=True
            )

    @requires_show
    def overlay(self, update_viewer: bool = True) -> None:
        """
        Re-index all poses or pdbstrings into a single `list` object at index `0`.

        Args:
            update_viewer: a `bool` object. If `True`, then update the visualization, which
                will automatically call the `show` method if it has not been called.
                If `False`, then do not update the visualization.
                Default: `True`

        Returns:
            `None`
        """
        for index in sorted(self.poses.keys()):
            if index != 0:
                self.poses[0].extend(self.poses[index])
                self.pdbstrings[0].extend(self.pdbstrings[index])
                self.poses.pop(index)
                self.pdbstrings.pop(index)
        if update_viewer:
            self.update_viewer(index=0)


def apply_metric_to_poses(
    metric: "PerResidueRealMetric", poses: Union[Iterable[Pose], Pose]
) -> None:
    _msg = "The 'poses' argument parameter must be a `Pose` object or an iterable of `Pose` objects. "
    if isinstance(poses, collections.abc.Iterable) and not isinstance(poses, Pose):
        for pose in poses:
            if isinstance(pose, Pose):
                with out:
                    metric.apply(pose)
            else:
                raise ValueError(_msg + f"Received: {type(pose)}")
    elif isinstance(poses, Pose):
        with out:
            metric.apply(poses)
    else:
        raise ValueError(_msg + f"Received: {type(poses)}")
