# coding: utf-8

from __future__ import absolute_import
from pyhml.models.typing import Typing
from .base_model_ import Model
from datetime import date, datetime
from typing import List, Dict
from ..util import deserialize_model


class Sample(Model):
    """

    Examples:

        >>> from pyhml.models.typing import Typing
        >>> from pyhml.models.sample import Sample

    """
    def __init__(self, center_code: int=None, id: str=None, collection_method: str=None, typing: List[Typing]=None):
        """
        Sample - a model defined in Swagger

        :param center_code: The center_code of this Sample.
        :type center_code: int
        :param id: The id of this Sample.
        :type id: str
        :param collection_method: The collection_method of this Sample.
        :type collection_method: str
        :param typing: The typing of this Sample.
        :type typing: List[Typing]
        """
        self.swagger_types = {
            'center_code': int,
            'id': str,
            'collection_method': str,
            'typing': List[Typing],
            'seq_records': Dict
        }

        self.attribute_map = {
            'center_code': 'center_code',
            'id': 'id',
            'collection_method': 'collection_method',
            'typing': 'typing',
            'seq_records': 'seq_records'
        }

        self._center_code = center_code
        self._id = id
        self._collection_method = collection_method
        self._typing = typing
        self._seq_records = {}

    @classmethod
    def from_dict(cls, dikt) -> 'Sample':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The Sample of this Sample.
        :rtype: Sample
        """
        return deserialize_model(dikt, cls)

    @property
    def center_code(self) -> int:
        """
        Gets the center_code of this Sample.

        :return: The center_code of this Sample.
        :rtype: int
        """
        return self._center_code

    @center_code.setter
    def center_code(self, center_code: int):
        """
        Sets the center_code of this Sample.

        :param center_code: The center_code of this Sample.
        :type center_code: int
        """

        self._center_code = center_code

    @property
    def id(self) -> str:
        """
        Gets the id of this Sample.

        :return: The id of this Sample.
        :rtype: str
        """
        return self._id

    @id.setter
    def id(self, id: str):
        """
        Sets the id of this Sample.

        :param id: The id of this Sample.
        :type id: str
        """

        self._id = id

    @property
    def collection_method(self) -> str:
        """
        Gets the collection_method of this Sample.

        :return: The collection_method of this Sample.
        :rtype: str
        """
        return self._collection_method

    @collection_method.setter
    def collection_method(self, collection_method: str):
        """
        Sets the collection_method of this Sample.

        :param collection_method: The collection_method of this Sample.
        :type collection_method: str
        """

        self._collection_method = collection_method

    @property
    def typing(self) -> List[Typing]:
        """
        Gets the typing of this Sample.

        :return: The typing of this Sample.
        :rtype: List[Typing]
        """
        return self._typing

    @typing.setter
    def typing(self, typing: List[Typing]):
        """
        Sets the typing of this Sample.

        :param typing: The typing of this Sample.
        :type typing: List[Typing]
        """

        self._typing = typing

    @property
    def seq_records(self) -> Dict:
        """
        Gets the seq_records of this Sample.

        :return: The seq_records of this Sample.
        :rtype: Dict
        """
        return self._seq_records

    def create_seqrecords(self):
        """
        Creates the seq_records for this Sample.

        :type seq_records: Dict
        """
        records = {}
        subid = self.id
        for typing in self.typing:
            loc, seq_rec = typing.create_seqrecord(subid)
            if loc not in records:
                records.update({loc: seq_rec})

        self._seq_records = records




