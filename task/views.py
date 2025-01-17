import json
import os
import pickle
import random
import time
import traceback
from pathlib import Path

from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.http import FileResponse
from rest_framework import viewsets
from rest_framework.decorators import api_view
from rest_framework.response import Response

import utils.analysis
from cell_type_workspace_api import settings_local as local_settings
from task.models import Task
from task.serializers import taskSerializer


class taskViewSet(viewsets.ModelViewSet):
    queryset = Task.objects.order_by('id')
    serializer_class = taskSerializer


@api_view(['GET'])
def task_detail_view(request):
    userid = request.query_params.dict()['userid']
    taskslist = Task.objects.filter(user=userid).order_by("-created_at")
    serializer = taskSerializer(taskslist, many=True)
    return Response({'results': serializer.data})


@api_view(['POST'])
def createtask(request):
    """
    Create a new task
    - userid
    - submitfile
    - taskname
    - tasktype
    - projectname
    - modulename
    - parameters
    """
    # create user task folder and save the file
    user_task_dir = request.data.get('taskdir') or str(int(time.time())) + '_' + str(random.randint(1000, 9999))
    user_path = local_settings.USER_TASK_PATH / user_task_dir
    upload_file_path = user_path / 'upload/'
    if request.data.get('taskdir'):
        if not os.path.exists(user_path):
            raise ValueError('taskdir not exists')
    else:
        os.makedirs(upload_file_path, exist_ok=False)
    file = request.FILES.get('submitfile')
    if file:
        default_storage.save((upload_file_path / 'input.zip').relative_to(local_settings.USER_TASK_PATH),
                             ContentFile(file.read()))

    # get parameters from request
    parameters_string = request.data['parameters']
    print(parameters_string)
    parameters_dict = json.loads(parameters_string)

    # create task object
    res = {}
    new_task = Task.objects.create(
        name=request.data['modulename'], user=request.data['userid'], userpath=user_task_dir,
        task_type=request.data['tasktype'], status='Created', modulelist=request.data['modulename'])

    # create module object and run the task
    if new_task.task_type == 'module':
        try:
            # run the task script
            def get_class_from_module(module, class_name):
                # 使用 getattr() 尝试从模块中获取类对象,如果类不存在，则返回 None
                return getattr(module, class_name, None)

            cls = get_class_from_module(utils.analysis, request.data['modulename'])

            if cls is None:
                res['status'] = 'Failed'
                new_task.status = 'Failed'
                res['message'] = 'module not found'
                raise ValueError('module not found')

            else:
                new_module = cls(request.data['taskname'], user_task_dir, parameters_dict)
                job_id = new_module.process()

                task_detail_json = [
                    {'modulename': request.data['modulename'], 'parameters_dict': parameters_dict, 'job_id': job_id,
                     'status': 'Created'}]
                with open(user_path / 'taskdetail.json', 'w') as f:
                    json.dump(task_detail_json, f, ensure_ascii=False, indent=4)
                with open(user_path / 'moduleobject.pkl', 'wb') as f:
                    pickle.dump(new_module, f)
                new_task.status = 'Running'
                res['status'] = 'Success'
                res['message'] = 'task create successfully'
                res['data'] = {'taskid': new_task.id}
        except Exception as e:
            res['status'] = 'Failed'
            res['message'] = e
            new_task.status = 'Failed'
            traceback.print_exc()
    new_task.save()
    return Response(res)


@api_view(['GET'])
def task_list_view(request):
    userid = request.query_params.dict()['userid']
    tasks_list = Task.objects.filter(user=userid).order_by("-created_at")
    serializer = taskSerializer(tasks_list, many=True)
    return Response({'results': serializer.data})


@api_view(['GET'])
def task_detail_view(request):
    task_id = request.query_params.dict()['taskid']
    task_object = Task.objects.filter(id=task_id)
    serializer = taskSerializer(task_object, many=True)
    task_data = serializer.data[0]
    task_data['userpath'] = task_data['userpath']
    task_data['inputpath'] = local_settings.FILEAPI + task_data['userpath'] + '/upload/input.csv'
    return Response({'results': task_data})


@api_view(['GET'])
def getoutputfile(request, path):
    file_path = local_settings.USER_TASK_PATH / path
    file = open(file_path, 'rb')
    response = FileResponse(file)
    filename = file.name.split('/')[-1]
    response['Content-Disposition'] = "attachment; filename=" + filename
    response['Content-Type'] = 'text/plain'
    return response


@api_view(['GET'])
def taskresultview(request):
    """
    Get task result: metadata, expression, umap, casuality
    - taskid*
    - resulttype*: metadata/expression/umap/batcheffect/casuality
    - metadata: page, pagesize (optional)
    - batcheffect: gene, compid (optional)
    - casuality: cluster (optional)
    """
    query_params = request.query_params.dict()
    taskid = query_params['taskid']
    taskobject = Task.objects.get(id=taskid)
    objectpath = local_settings.USER_TASK_PATH / taskobject.userpath / 'moduleobject.pkl'
    with open(objectpath, 'rb') as f:
        # 载入模块对象
        module = pickle.load(f)
    res = module.getresult(query_params)
    return Response(res)
