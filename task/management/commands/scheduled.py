import datetime
import json
import pickle

from django.core.management.base import BaseCommand, CommandError

from cell_type_workspace_api import settings_local as local_settings
from cell_type_workspace_api.settings import BASE_DIR
from task.models import Task


class Command(BaseCommand):
    help = 'Describe the command here'

    def handle(self, *args, **options):
        self.stdout.write('Starting script...')
        try:
            tasklist = Task.objects.filter(status='Running')
            for task in tasklist:
                if task.task_type == 'module':
                    objectpath = local_settings.USER_TASK_PATH / task.userpath / 'moduleobject.pkl'
                    with open(objectpath, 'rb') as f:
                        taskobject = pickle.load(f)

                    if taskobject.check_status() != 'COMPLETED':
                        continue
                    else:
                        task.status = 'Completed'
                        task.save()
                        with open(local_settings.USER_TASK_PATH / task.userpath / 'taskdetail.json', 'r') as f:
                            jsondata = json.load(f)
                        jsondata[0]['status'] = 'Completed'
                        with open(local_settings.USER_TASK_PATH / task.userpath / 'taskdetail.json', 'w') as f:
                            json.dump(jsondata, f, ensure_ascii=False, indent=4)
                        with open(local_settings.USER_TASK_PATH / task.userpath / 'moduleobject.pkl', 'wb') as f:
                            pickle.dump(taskobject, f)
            current_time = datetime.datetime.now()
            with open(BASE_DIR / 'workspace/log/update.txt', 'a+') as f:
                f.write('exec update start  ' + str(current_time) + "\n")
        except Exception as e:
            raise CommandError('Something went wrong: %s' % e)
            self.stdout.write(self.style.ERROR('Error occurred'))
